import argparse
from pathlib import Path
import math
import uuid
import copy
from datetime import datetime, date
import shutil
import os
import re

# CAD libraries
import ezdxf
import numpy as np
from shapely.geometry import Polygon, LineString, Point as ShapelyPoint, MultiPolygon, MultiPoint
from shapely.ops import unary_union, polygonize

# GIS & Data libraries
try:
    import geopandas as gpd
    import fiona
    import pandas as pd
except ImportError:
    print("LET OP: 'geopandas', 'fiona' of 'pandas' is niet geïnstalleerd.")
    gpd = None

# ==============================================================================
# CONFIGURATIE & OTL MAPPING LOGICA
# ==============================================================================
OTL_LAYER_MAPPING = {
    "BSS": "Betonsteen",
    "TEGEL": "Betontegel",
    "GEBAKKEN": "Gebakken steen",
    "SBS": "Gebakken steen",
    "GRIND": "Grind",
    "GRINDTEGEL": "Grindtegel",
    "NATUURSTEEN": "Natuursteentegel",      
    "GRAVEL": "Gravel", 
    "PUIN": "Puin",
    "BASALT": "Basaltblokken", 
    "BETONPLAAT": "Betonplaat", 
    "KEI": "Kei",
    "ASFALT": "Asfaltbeton",    
    "ZAND": "Zand",
    "GRAS": "Gras verhard",
    "HOUT": "Houtsnippers"
}
DEFAULT_LAYER = "Open grond"

POINT_LAYER_MAPPING = {
    "TROTTOIR": "Trottoirkolk",
    "POMPPUT": "Opvangput",
    "STRAAT": "Straatkolk",
    "GOOT": "Roostergoot",
    "INFILTRATIE": "Infiltratiekolk"
}

SYSTEM_MAPPING = {
    "identificatie": "Identificatie",       
    "imgeo_id": "IMGeo-identificatie", 
    "jaar_aanleg": "Jaar van aanleg",
    "jaar_uitgevoerd_onderhoud": "Jaar uitgevoerd onderhoud",  
    "opleverdatum": "Opleverdatum",  
    "begin_garantie": "Begin garantieperiode",
    "einde_garantie": "Einde garantieperiode",
    "bronhouder": "Bronhouder",              
    "materiaal": "Materiaal",
    "kleur": "Kleur",
    "functie": "Gebruiksfunctie"
}

COLUMN_MAPPING = {
    "MATERIAAL": "Materiaal",
    "TYPE": "Type rijstrook", 
    "KLEUR": "Kleur",
    "STATUS": "Status",
    "FUNCTIE": "Gebruiksfunctie",
    "AFMETING": "Afmeting",
    "FORMAAT": "Formaat",
    "FABRIKANT": "Fabrikant"
}

VALUE_MAPPING = {
    "RD": "Rood",
    "RB": "Rood bruin",
    "B": "Bruin",
    "R": "Rood",
    "G": "Grijs",
    "ZW": "Zwart",
    "KF": "Keiformaat",
}

# --- NIEUWE LOGICA VOOR CONSTRUCTIELAGEN ---
MATERIAL_CATEGORIES = {
    "Asfaltbeton": "ASFALT",
    "Steenmastiekasfalt": "ASFALT",
    "Printasfalt": "ASFALT",
    "Asfaltbeton dunne laag": "ASFALT",
    "Gebakken steen": "ELEMENTEN",
    "Betonsteen": "ELEMENTEN",
    "Natuursteenkei": "ELEMENTEN",
    "Kei": "ELEMENTEN",
    "Basaltblokken": "ELEMENTEN",
    "Betontegel": "TEGELS",
    "Natuursteentegel": "TEGELS",
    "Grindtegel": "TEGELS",
    "Rubbertegel": "TEGELS",
    "Grind": "HALFVERHARDING",
    "Gravel": "HALFVERHARDING",
    "Schelpen": "HALFVERHARDING",
    "Puin": "HALFVERHARDING",
    "Gralux": "HALFVERHARDING",
    "Houtsnippers": "HALFVERHARDING",
    "Boomschors": "HALFVERHARDING",
    "Betonplaat": "BETON",
    "Gewapend beton": "BETON",
    "Verdeuveld": "BETON",
    "Onverdeuveld": "BETON"
}

# Opbouw per categorie
LAYER_BUILDUP_RULES = {
    "ASFALT": [
        (1, "DEKLAAG", 40, "Asfalt deklaag"),
        (2, "TUSSENLAAG", 60, "AC binder"),
        (3, "ONDERLAAG", 100, "AC base"),
        (4, "FUNDERING", 250, "Menggranulaat"),
        (5, "ZANDBED", 200, "Zand")
    ],
    "ELEMENTEN": [
        (1, "ELEMENTENLAAG", 80, "Bestrating"),
        (2, "LEGBED", 40, "Straatzand"),
        (3, "FUNDERING", 250, "Menggranulaat"),
        (4, "ZANDBED", 200, "Zand")
    ],
    "TEGELS": [
        (1, "ELEMENTENLAAG", 50, "Tegels"),
        (2, "LEGBED", 40, "Straatzand"),
        (3, "ZANDBED", 200, "Zand")
    ],
    "HALFVERHARDING": [
        (1, "DEKLAAG", 100, "Halfverharding"),
        (2, "ZANDBED", 200, "Zand")
    ],
    "BETON": [
        (1, "DEKLAAG", 200, "Beton"),
        (2, "FUNDERING", 200, "Menggranulaat"),
        (3, "ZANDBED", 200, "Zand")
    ],
    "DEFAULT": [
        (1, "DEKLAAG", 50, "Onbekend"),
        (2, "ZANDBED", 200, "Zand")
    ]
}

# ==============================================================================
# GEOMETRIE & CAD LOGICA
# ==============================================================================
ARC_MAX_DEG = 5.0

def snap_endpoints(lines, layers, tolerance):
    all_endpoints = [] 
    endpoint_to_line = []
    for line_idx, line in enumerate(lines):
        coords = list(line.coords)
        layer = layers[line_idx]
        all_endpoints.append((coords[0][0], coords[0][1], layer))
        endpoint_to_line.append((line_idx, "start"))
        all_endpoints.append((coords[-1][0], coords[-1][1], layer))
        endpoint_to_line.append((line_idx, "end"))
    point_mapping = {}
    used_points = set()
    for i, (x1, y1, layer1) in enumerate(all_endpoints):
        if i in used_points: continue
        cluster = [(x1, y1)]
        cluster_indices = [i]
        used_points.add(i)
        for j, (x2, y2, layer2) in enumerate(all_endpoints):
            if j in used_points: continue
            if layer2 != layer1: continue
            dist = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
            if dist < tolerance:
                cluster.append((x2, y2))
                cluster_indices.append(j)
                used_points.add(j)
        avg_x = sum(p[0] for p in cluster) / len(cluster)
        avg_y = sum(p[1] for p in cluster) / len(cluster)
        snapped_pos = (avg_x, avg_y)
        for idx in cluster_indices:
            point_mapping[idx] = snapped_pos
    snapped_lines = []
    for line_idx, line in enumerate(lines):
        coords = list(line.coords)
        start_endpoint_idx = line_idx * 2
        end_endpoint_idx = line_idx * 2 + 1
        start_snapped = point_mapping.get(start_endpoint_idx, coords[0])
        end_snapped = point_mapping.get(end_endpoint_idx, coords[-1])
        new_coords = [start_snapped] + coords[1:-1] + [end_snapped]
        snapped_lines.append(LineString(new_coords))
    return snapped_lines

def apply_cutout_holes(polygons, blocks_data, min_hole_area=0.05):
    def poly_has_block(poly):
        for b in blocks_data:
            if poly.contains(b["location"]): return True
        return False
    polygons_with_svh = []
    for poly in polygons:
        has_svh = False
        for block in blocks_data:
            if block["name"].startswith("SVH") and poly.contains(block["location"]):
                has_svh = True
                break
        if has_svh: polygons_with_svh.append(poly)
        else: polygons_with_svh.append(poly)
    sorted_polys = sorted(polygons_with_svh, key=lambda p: p.area, reverse=True)
    used_as_holes = set()
    result_polygons = []
    for i, large_poly in enumerate(sorted_polys):
        if i in used_as_holes: continue
        holes = []
        for j, small_poly in enumerate(sorted_polys):
            if i == j or j in used_as_holes: continue
            if not large_poly.contains(small_poly): continue
            if small_poly.area < min_hole_area and not poly_has_block(small_poly): continue
            holes.append(small_poly)
            used_as_holes.add(j)
        if holes:
            hole_coords = [list(hole.exterior.coords) for hole in holes]
            new_poly = Polygon(list(large_poly.exterior.coords), holes=hole_coords)
            result_polygons.append(new_poly)
        else:
            result_polygons.append(large_poly)
    return result_polygons

def polyline_to_linestring(ent, arc_max_deg=ARC_MAX_DEG):
    dxftype = ent.dxftype()
    has_arc = bool(getattr(ent, "has_arc", False))
    if not has_arc:
        if dxftype == "LWPOLYLINE":
            pts = [(p[0], p[1]) for p in ent.get_points()]
            if getattr(ent, "closed", False):
                if pts and pts[0] != pts[-1]: pts.append(pts[0])
            if len(pts) >= 2: return LineString(pts)
            return None
        if dxftype == "POLYLINE":
            pts = []
            for v in ent.vertices:
                if hasattr(v.dxf, "location"):
                    pts.append((v.dxf.location.x, v.dxf.location.y))
            if getattr(ent, "is_closed", False):
                if pts and pts[0] != pts[-1]: pts.append(pts[0])
            if len(pts) >= 2: return LineString(pts)
            return None
    pts = []
    try: segments = list(ent.virtual_entities())
    except Exception: segments = []
    for seg in segments:
        t = seg.dxftype()
        if t == "LINE":
            start = seg.dxf.start
            end = seg.dxf.end
            if not pts: pts.append((start.x, start.y))
            pts.append((end.x, end.y))
        elif t == "ARC":
            center = seg.dxf.center
            r = seg.dxf.radius
            start_angle = math.radians(seg.dxf.start_angle)
            end_angle = math.radians(seg.dxf.end_angle)
            total_angle = end_angle - start_angle
            if total_angle > math.pi: total_angle -= 2 * math.pi
            elif total_angle < -math.pi: total_angle += 2 * math.pi
            if total_angle == 0: continue
            max_step = math.radians(arc_max_deg)
            n = max(4, int(abs(total_angle) / max_step))
            for i in range(n + 1):
                a = start_angle + total_angle * (i / n)
                x = center.x + r * math.cos(a)
                y = center.y + r * math.sin(a)
                if not pts or (x, y) != pts[-1]: pts.append((x, y))
    if getattr(ent, "closed", False) or getattr(ent, "is_closed", False):
        if pts and pts[0] != pts[-1]: pts.append(pts[0])
    if len(pts) < 2: return None
    return LineString(pts)

def read_lines_and_blocks(doc):
    msp = doc.modelspace()
    lines_with_props = []
    for line in msp.query("LINE"):
        try:
            start = line.dxf.start
            end = line.dxf.end
            geom = LineString([(start.x, start.y), (end.x, end.y)])
            lines_with_props.append({"geometry": geom, "entity": line, "type": "LINE"})
        except Exception: pass
    for pl in msp.query("LWPOLYLINE"):
        try:
            geom = polyline_to_linestring(pl)
            if geom is not None:
                lines_with_props.append({"geometry": geom, "entity": pl, "type": "LWPOLYLINE"})
        except Exception: pass
    for pl in msp.query("POLYLINE"):
        try:
            geom = polyline_to_linestring(pl)
            if geom is not None:
                lines_with_props.append({"geometry": geom, "entity": pl, "type": "POLYLINE"})
        except Exception: pass
    blocks_data = []
    for insert in msp.query("INSERT"):
        block_name = insert.dxf.name if hasattr(insert.dxf, "name") else str(insert)
        try:
            location = insert.dxf.insert
            if location is None: raise AttributeError
            shapely_location = ShapelyPoint(location.x, location.y)
        except Exception: continue
        block_color = None
        if hasattr(insert.dxf, "color") and insert.dxf.color != 256:
            block_color = insert.dxf.color
        if block_color is None:
             try:
                layer_name = insert.dxf.layer
                if layer_name in doc.layers:
                    layer = doc.layers.get(layer_name)
                    if hasattr(layer.dxf, "color"): block_color = layer.dxf.color
             except: pass
        attributes = {}
        try:
            for attrib in insert.attribs:
                tag = attrib.dxf.tag
                text = attrib.dxf.text
                if tag and text: attributes[tag] = text
        except Exception as e: pass
        blocks_data.append(
            {
                "name": block_name, 
                "location": shapely_location, 
                "color": block_color, 
                "layer": insert.dxf.layer, 
                "insert": insert, 
                "attributes": attributes
            }
        )
    return lines_with_props, blocks_data

def build_polygons(lines_with_props, tolerance):
    lines_geom = [item["geometry"] for item in lines_with_props]
    if not lines_geom: return []
    layers = [item["entity"].dxf.layer for item in lines_with_props]
    snapped = snap_endpoints(lines_geom, layers, tolerance)
    merged = unary_union(snapped)
    polys = list(polygonize(merged))
    polygons = [p for p in polys if p.is_valid and p.area > 0.01]
    return polygons

def classify_polygons(polygons, blocks_data):
    svh_polygons = []
    non_svh_polygons = []
    for poly in polygons:
        material = None
        block_info = None
        for block in blocks_data:
            if poly.contains(block["location"]):
                if block["name"].startswith("SVH"):
                    material = block["name"]
                    block_info = block
                    break
        if material:
            svh_polygons.append({"polygon": poly, "material": material, "area": poly.area, "block_info": block_info})
        else:
            non_svh_polygons.append({"polygon": poly, "area": poly.area})
    return svh_polygons, non_svh_polygons

# ==============================================================================
# GIS EXPORT LOGICA
# ==============================================================================

def prepare_full_gpkg_copy(input_dxf_path, template_path):
    if not os.path.exists(template_path):
        print(f"FOUT: Template niet gevonden: {template_path}")
        return None
    input_path = Path(input_dxf_path)
    stem = input_path.stem 
    project_name = stem.split('_')[0]
    new_filename = f"{project_name}_GEODATA_d5.0.gpkg"
    gpkg_path = input_path.parent / new_filename
    
    if gpkg_path.exists():
        try:
            os.remove(gpkg_path)
        except Exception as e:
            print(f"Kon oud bestand niet verwijderen: {e}")
            return None
    try:
        shutil.copy(template_path, gpkg_path)
    except Exception as e:
        print(f"FOUT bij kopiëren template: {e}")
        return None
    # Verwijder eventuele administratieve 'Metadata' laag uit de kopie
    try:
        existing = fiona.listlayers(str(gpkg_path))
        for lyr in existing:
            if lyr.lower() == 'metadata':
                try:
                    fiona.remove(str(gpkg_path), layer=lyr)
                    print(f"Verwijderd administratieve laag uit export: {lyr}")
                except Exception as e:
                    print(f"Kon laag '{lyr}' niet verwijderen uit kopie: {e}")
    except Exception:
        pass
    return str(gpkg_path)

def clean_up_gpkg(gpkg_path):
    if not os.path.exists(gpkg_path): return
    layers_to_remove = []
    try:
        all_layers = fiona.listlayers(gpkg_path)
        for layer_name in all_layers:
            # Sla de administratieve tabellen over met schoonmaken!
            if layer_name.lower() in ["constructielaag", "metadata"]: continue
            try:
                with fiona.open(gpkg_path, layer=layer_name) as src:
                    if len(src) == 0:
                        layers_to_remove.append(layer_name)
            except: pass
    except Exception as e:
        print(f"Fout bij scannen lagen: {e}")
        return
        
    if layers_to_remove:
        print(f"--- Opschonen: {len(layers_to_remove)} lege lagen verwijderen ---")
        for layer in layers_to_remove:
            try:
                fiona.remove(gpkg_path, layer=layer)
            except Exception as e:
                print(f"Kon laag '{layer}' niet verwijderen: {e}")

def get_mapped_attributes(cad_attributes):
    mapped = {}
    for tag, otl_col in COLUMN_MAPPING.items():
        if tag in cad_attributes:
            raw_val = cad_attributes[tag]
            final_val = VALUE_MAPPING.get(raw_val.upper(), raw_val)
            mapped[otl_col] = final_val
    return mapped

def write_to_layer(gpkg_path, layer_name, records):
    """
    Robsuste schrijffunctie die Schema's uitlijnt met de template
    zodat GeoPandas niet per ongeluk kolommen op de verkeerde plek wegschrijft.
    """
    if not records: return
    print(f" -> Verwerken laag: {layer_name} ({len(records)} objecten)...")
    
    has_geom = any(r.get('geometry') is not None for r in records)
    
    if has_geom:
        from shapely.geometry import shape
        normalized = []
        for r in records:
            geom = r.get('geometry')
            if geom is None:
                normalized.append(r)
                continue
            # If geometry is a geo-interface (dict), convert to shapely
            try:
                if not hasattr(geom, 'geom_type'):
                    geom = shape(geom)
            except Exception:
                pass

            gtype = getattr(geom, 'geom_type', '').lower()
            if gtype == 'polygon':
                try:
                    geom = MultiPolygon([geom])
                except Exception:
                    pass
            elif gtype == 'point':
                try:
                    geom = MultiPoint([geom])
                except Exception:
                    pass

            r['geometry'] = geom
            normalized.append(r)

        gdf = gpd.GeoDataFrame(normalized, crs="EPSG:28992")
    else:
        # Als we geen geometry hebben (bijv Constructielaag)
        gdf = pd.DataFrame(records)
        gdf['geometry'] = None
        gdf = gpd.GeoDataFrame(gdf, geometry='geometry', crs="EPSG:28992")

    # Opschonen van haakjes in ALLE string velden
    for col in gdf.columns:
        if col != 'geometry' and gdf[col].dtype == object:
            gdf[col] = gdf[col].apply(lambda x: re.sub(r'\s*\([^)]*\)', '', str(x)).strip() if isinstance(x, str) and pd.notnull(x) else x)

    # SCHEMA UITLIJNING: Fix de "2026 in Bijzondere constructie" warning.
    try:
        with fiona.open(gpkg_path, layer=layer_name) as target_info:
            target_cols = list(target_info.schema['properties'].keys())
        
        # 1. Voeg ontbrekende target kolommen leeg toe aan GDF
        for c in target_cols:
            if c not in gdf.columns:
                gdf[c] = None
        
        # 2. Pak exact de kolommen uit de target en in die exacte volgorde
        cols_to_keep = [c for c in target_cols if c in gdf.columns]
        if 'geometry' in gdf.columns:
            cols_to_keep.append('geometry')
            
        gdf = gdf[cols_to_keep]
    except Exception as e:
        print(f"    Schema alignment waarschuwing: {e}")

    # Exporteer
    try:
        gdf.to_file(gpkg_path, layer=layer_name, driver="GPKG", mode="a")
        print(f"    Succes via GeoPandas!")
    except Exception as e:
        print(f"    GeoPandas opslag faalde ({e}). Fallback naar Fiona wordt gestart...")
        try:
            with fiona.open(gpkg_path, mode='a', layer=layer_name) as dst:
                schema_cols = list(dst.schema['properties'].keys())
                schema_map = {col.lower(): col for col in schema_cols}
                fiona_records = []
                
                for i, row in gdf.iterrows():
                    geom_val = None
                    if has_geom and pd.notnull(row.get('geometry')):
                        geom_val = row['geometry'].__geo_interface__
                    
                    props = {}
                    for script_col, value in row.items():
                        if script_col == 'geometry': continue
                        target_col = schema_map.get(script_col.lower())
                        if target_col:
                            props[target_col] = None if pd.isna(value) else value
                            
                    for real_col in schema_cols:
                        if real_col not in props:
                            props[real_col] = None
                    
                    fiona_records.append({'geometry': geom_val, 'properties': props})
                
                dst.writerecords(fiona_records)
                print(f"    Succes via Fiona fallback!")
        except Exception as e2:
            print(f"    FOUT bij schrijven naar {layer_name} (Fiona faalde ook): {e2}")

def export_svh_to_gis(svh_polygons, gpkg_path, project_config):
    if not svh_polygons or not gpkg_path: return
    layers_data = {}
    constructielagen = [] 
    
    for item in svh_polygons:
        poly = item["polygon"]
        block_info = item["block_info"]
        attrs = block_info.get("attributes", {})
        material_name = item["material"].upper()
        
        target_layer = DEFAULT_LAYER
        for keyword, layer_name in OTL_LAYER_MAPPING.items():
            if keyword in material_name:
                target_layer = layer_name
                break
        
        poly_identificatie = str(uuid.uuid4())
        
        internal_record = {
            "geometry": poly,
            "identificatie": poly_identificatie, 
            "imgeo_id": str(uuid.uuid4()),      
            "jaar_aanleg": project_config.get("jaar_aanleg"),
            "jaar_uitgevoerd_onderhoud": project_config.get("jaar_uitgevoerd_onderhoud"),
            "opleverdatum": project_config.get("opleverdatum"),
            "begin_garantie": project_config.get("begin_garantie"),
            "einde_garantie": project_config.get("einde_garantie"),
            "bronhouder": project_config.get("bronhouder"),
            "materiaal": item["material"] 
        }
        mapped_attrs = get_mapped_attributes(attrs)
        internal_record.update(mapped_attrs)
        
        final_record = {"geometry": poly}
        for k, v in internal_record.items():
            if k == "geometry": continue
            target_key = SYSTEM_MAPPING.get(k, k) 
            final_record[target_key] = v

        if target_layer not in layers_data: layers_data[target_layer] = []
        layers_data[target_layer].append(final_record)

        # --- AANMAKEN CONSTRUCTIELAGEN O.B.V. TYPE ---
        categorie = MATERIAL_CATEGORIES.get(target_layer, "DEFAULT")
        opbouw_regels = LAYER_BUILDUP_RULES.get(categorie, LAYER_BUILDUP_RULES["DEFAULT"])
        
        for (laag_nr, soort, default_dikte, sub_materiaal) in opbouw_regels:
            if laag_nr == 1:
                dikte_str = attrs.get("DIKTE", str(default_dikte))
                try: dikte = int(float(dikte_str))
                except: dikte = default_dikte
                mat_naam = mapped_attrs.get("Materiaal", item["material"])
            else:
                dikte = default_dikte
                mat_naam = sub_materiaal
                
            constructielaag_record = {
                "geometry": None, 
                "Identificatie": str(uuid.uuid4()),
                "ID-Verhardingsobject": poly_identificatie, 
                "Constructielaagsoort": soort,
                "Dikte constructielaag": dikte,
                "Jaar van aanleg": project_config.get("jaar_aanleg"),
                "Laagnummer": laag_nr,
                "Materiaal": mat_naam
            }
            constructielagen.append(constructielaag_record)

    print("--- Start export VLAKKEN ---")
    for layer_name, records in layers_data.items():
        write_to_layer(gpkg_path, layer_name, records)

    if constructielagen:
        print("--- Start export CONSTRUCTIELAGEN ---")
        write_to_layer(gpkg_path, "Constructielaag", constructielagen)

def export_sri_to_gis(sri_blocks, gpkg_path, project_config):
    if not sri_blocks or not gpkg_path: return
    layers_data = {}
    
    for block in sri_blocks:
        point = block["location"] 
        attrs = block.get("attributes", {})
        block_name = block["name"].upper()
        
        target_layer = None
        for keyword, layer_name in POINT_LAYER_MAPPING.items():
            if keyword in block_name:
                target_layer = layer_name
                break
        if target_layer is None: continue
        
        internal_record = {
            "geometry": point,
            "identificatie": str(uuid.uuid4()),
            "jaar_aanleg": project_config.get("jaar_aanleg"),
            "jaar_uitgevoerd_onderhoud": project_config.get("jaar_uitgevoerd_onderhoud"),
            "opleverdatum": project_config.get("opleverdatum"),
            "begin_garantie": project_config.get("begin_garantie"),
            "einde_garantie": project_config.get("einde_garantie"),
            "bronhouder": project_config.get("bronhouder"),
        }
        mapped_attrs = get_mapped_attributes(attrs)
        internal_record.update(mapped_attrs)
        
        final_record = {"geometry": point}
        for k, v in internal_record.items():
            if k == "geometry": continue
            target_key = SYSTEM_MAPPING.get(k, k)
            final_record[target_key] = v
        
        if target_layer not in layers_data: layers_data[target_layer] = []
        layers_data[target_layer].append(final_record)

    print("--- Start export PUNTEN ---")
    for layer_name, records in layers_data.items():
        write_to_layer(gpkg_path, layer_name, records)

# ==============================================================================
# HOEVEELHEDEN EXPORT (EXCEL)
# ==============================================================================
def export_quantities_to_excel(output_dxf_path, project_config, svh_polygons, lines_with_props, blocks_data):
    try:
        input_path = Path(output_dxf_path)
        stem = input_path.stem.replace("_verwerkt", "") 
        project_name = stem.split('_')[0]
        excel_filename = f"{project_name}_Hoeveelheden.xlsx"
        excel_path = input_path.parent / excel_filename

        data_rows = []

        data_rows.append({"Categorie": "PROJECT INFORMATIE", "Omschrijving": "", "Hoeveelheid": "", "Eenheid": ""})
        data_rows.append({"Categorie": "", "Omschrijving": "Project", "Hoeveelheid": project_name, "Eenheid": ""})
        data_rows.append({"Categorie": "", "Omschrijving": "Datum", "Hoeveelheid": datetime.now().strftime("%d-%m-%Y"), "Eenheid": ""})
        data_rows.append({"Categorie": "", "Omschrijving": "", "Hoeveelheid": "", "Eenheid": ""}) 

        data_rows.append({"Categorie": "VERHARDING (VLAKKEN)", "Omschrijving": "", "Hoeveelheid": "", "Eenheid": ""})
        
        area_stats = {}
        for poly_data in svh_polygons:
            # RegEx filter toepassen in Excel output voor strakke presentatie
            raw_mat = poly_data["material"]
            mat = re.sub(r'\s*\([^)]*\)', '', raw_mat).strip()
            
            area_stats[mat] = area_stats.get(mat, 0.0) + poly_data["area"]
        
        for mat, area in area_stats.items():
            data_rows.append({
                "Categorie": "Verharding",
                "Omschrijving": mat,
                "Hoeveelheid": round(area, 2),
                "Eenheid": "m2"
            })

        data_rows.append({"Categorie": "", "Omschrijving": "", "Hoeveelheid": "", "Eenheid": ""}) 

        data_rows.append({"Categorie": "INRICHTING (PUNTEN)", "Omschrijving": "", "Hoeveelheid": "", "Eenheid": ""})
        
        block_stats = {}
        for block in blocks_data:
            name = block["name"]
            layer = block.get("layer", "").upper()
            
            if name.startswith("SVH"): continue
            if layer.startswith("X"): continue
            
            # RegEx filter toepassen voor Inrichting
            clean_name = re.sub(r'\s*\([^)]*\)', '', name).strip()
            block_stats[clean_name] = block_stats.get(clean_name, 0) + 1
            
        for name, count in block_stats.items():
            data_rows.append({
                "Categorie": "Inrichting",
                "Omschrijving": name,
                "Hoeveelheid": count,
                "Eenheid": "st"
            })

        df = pd.DataFrame(data_rows)
        df.to_excel(excel_path, index=False, engine='xlsxwriter')
        
        print(f"Excel aangemaakt: {excel_path}")
        return str(excel_path)

    except Exception as e:
        print(f"Fout bij maken Excel: {e}")
        return None

# ==============================================================================
# DXF EXPORT
# ==============================================================================
def export_dxf(doc, svh_polygons, non_svh_polygons, output_path, color_palette=None):
    if not svh_polygons and not non_svh_polygons:
        raise RuntimeError("Geen vlakken gevonden om te exporteren.")
    if color_palette is None:
        color_palette = [2, 3, 4, 5, 6, 30, 40, 50, 140, 150, 160, 170, 180, 190, 200]
    new_doc = doc
    msp = new_doc.modelspace()
    APP_ID = "SVH_DATA"
    if not new_doc.appids.has_entry(APP_ID): new_doc.appids.add(APP_ID)
    material_colors = {}
    color_index = 0
    unique_materials = sorted(set(p["material"] for p in svh_polygons))
    for mat in unique_materials:
        material_colors[mat] = color_palette[color_index % len(color_palette)]
        color_index += 1
        layer_name = f"HATCH_{mat}".replace(" ", "_")
        if layer_name not in new_doc.layers:
            new_doc.layers.add(layer_name, color=material_colors[mat])
    if non_svh_polygons and "HATCH_GEEN_SVH" not in new_doc.layers:
        new_doc.layers.add("HATCH_GEEN_SVH", color=1)

    def create_hatch_entity(msp, poly, layer, block_data=None):
        coords = list(poly.exterior.coords)
        pts = [(x, y, 0) for x, y in coords]
        hatch = msp.add_hatch(color=ezdxf.const.BYLAYER, dxfattribs={"layer": layer})
        hatch.paths.add_polyline_path(pts, is_closed=True)
        for interior in poly.interiors:
            hole_pts = [(x, y, 0) for x, y in interior.coords]
            hatch.paths.add_polyline_path(hole_pts, is_closed=True)
        hatch.set_solid_fill()
        if block_data and "attributes" in block_data:
            xdata = []
            for tag, val in block_data["attributes"].items():
                xdata.append((1000, f"{tag}={val}"))
            xdata.append((1000, f"SOURCE_BLOCK={block_data['name']}"))
            if xdata: hatch.set_xdata(APP_ID, xdata)
        return hatch

    for poly_data in svh_polygons:
        poly = poly_data["polygon"]
        material = poly_data["material"]
        area = poly_data["area"]
        block_info = poly_data["block_info"]
        layer_name = f"HATCH_{material}".replace(" ", "_")
        try: create_hatch_entity(msp, poly, layer_name, block_info)
        except Exception: pass
        text_lines = [f"{material}", f"{area:.2f} m2"]
        if block_info and "attributes" in block_info:
            attrs = block_info["attributes"]
            priority_keys = ["FUNCTIE", "FORMAAT", "KLEUR", "TYPE"] 
            for key in priority_keys:
                if key in attrs: text_lines.append(f"{key}: {attrs[key]}")
        full_text = "\\P".join(text_lines)
        centroid = poly.centroid
        msp.add_mtext(full_text, dxfattribs={"layer": layer_name, "color": ezdxf.const.BYLAYER, "char_height": 0.5, "insert": (centroid.x, centroid.y, 0), "attachment_point": 5})

    for poly_data in non_svh_polygons:
        poly = poly_data["polygon"]
        area = poly_data["area"]
        layer_name = "HATCH_GEEN_SVH"
        try: create_hatch_entity(msp, poly, layer_name, None)
        except Exception: pass
        centroid = poly.centroid
        text = f"GEEN SVH\\P{area:.2f} m2"
        msp.add_mtext(text, dxfattribs={"layer": layer_name, "color": ezdxf.const.BYLAYER, "char_height": 0.5, "insert": (centroid.x, centroid.y, 0), "attachment_point": 5})

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    new_doc.saveas(str(output_path))
    return output_path

# ==============================================================================
# MAIN
# ==============================================================================
def process_dxf(input_path, output_path=None, tolerance=0.5, template_path="OTL_Verhardingen_v2_6_2.gpkg", user_metadata=None):
    input_path = Path(input_path)
    if output_path is None: output_path = input_path.with_name(input_path.stem + "_verwerkt.dxf")
    if not input_path.exists(): raise FileNotFoundError(f"Input niet gevonden: {input_path}")
    
    oplever_str = datetime.now().strftime("%Y-%m-%d")
    begin_garantie_str = oplever_str
    einde_garantie_str = oplever_str
    jaar_aanleg_val = str(datetime.now().year)

    if user_metadata:
        if user_metadata.get('opleverdatum'):
            dt_oplever = pd.to_datetime(user_metadata['opleverdatum'])
            dt_einde = dt_oplever + pd.DateOffset(months=6)
            oplever_str = dt_oplever.strftime("%Y-%m-%d")
            begin_garantie_str = oplever_str
            einde_garantie_str = dt_einde.strftime("%Y-%m-%d")
            jaar_aanleg_val = str(dt_oplever.year)

    project_config = {
        "bronhouder": "Gemeente Amsterdam",
        "jaar_aanleg": jaar_aanleg_val,
        "jaar_uitgevoerd_onderhoud": jaar_aanleg_val,
        "opleverdatum": oplever_str,
        "begin_garantie": begin_garantie_str,
        "einde_garantie": einde_garantie_str
    }

    print(f"Lezen DXF: {input_path}")
    doc = ezdxf.readfile(str(input_path))
    lines_with_props, blocks_data = read_lines_and_blocks(doc)
    
    sri_blocks = [b for b in blocks_data if "SRI" in b["name"]]
    
    polygons = build_polygons(lines_with_props, tolerance=tolerance)
    if not polygons: raise RuntimeError("Geen gesloten vlakken gevonden.")
    polygons_cut = apply_cutout_holes(polygons, blocks_data, min_hole_area=0.05)
    svh_polygons, non_svh_polygons = classify_polygons(polygons_cut, blocks_data)
    
    if gpd:
        gpkg_path = prepare_full_gpkg_copy(output_path, template_path)
        if gpkg_path:
            export_svh_to_gis(svh_polygons, gpkg_path, project_config)
            export_sri_to_gis(sri_blocks, gpkg_path, project_config)
            clean_up_gpkg(gpkg_path)
            
        export_quantities_to_excel(output_path, project_config, svh_polygons, lines_with_props, blocks_data)
            
    export_path = export_dxf(doc, svh_polygons, non_svh_polygons, output_path)
    print(f"\nKlaar. CAD output: {export_path}")

def main():
    parser = argparse.ArgumentParser(description="Revisie naar GIS (OTL Punten & Vlakken) + CAD.")
    parser.add_argument("input", help="Pad naar DXF")
    parser.add_argument("-o", "--output", help="Output pad", default=None)
    parser.add_argument("-t", "--tolerance", type=float, default=0.1) 
    parser.add_argument("--template", help="Pad naar OTL GPKG", default="OTL_Verhardingen_v2_6_2.gpkg")
    args = parser.parse_args()
    process_dxf(args.input, args.output, args.tolerance, args.template)

if __name__ == "__main__":
    main()
