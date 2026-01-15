import streamlit as st
import os
import shutil
import tempfile
import geopandas as gpd
import folium
from streamlit_folium import st_folium
import fiona
import ezdxf
from pyproj import Transformer
from backend_logic import process_dxf 
from PIL import Image
from datetime import date

# --- 1. CONFIGURATIE (MOET ALS EERSTE) ---
st.set_page_config(
    page_title="Kepa Revisie Portal", 
    page_icon="📐", 
    layout="wide", 
    initial_sidebar_state="expanded" 
)

# --- 2. CSS STYLING (HUISSTIJL & CORAL) ---
st.markdown("""
    <style>
        :root {
            --primary-color: #FF7F50; /* CORAL */
            --primary-hover: #E06A3D;
            --text-dark: #1d2b36;
            --bg-light: #f8f9fa;
        }

        #MainMenu {visibility: hidden;}
        footer {visibility: hidden;}
        
        [data-testid="stHeader"] {
            background-color: transparent;
        }
        [data-testid="stDecoration"] {
            display: none;
        }

        [data-testid="stSidebar"] {
            background-color: var(--bg-light);
            border-right: 1px solid #e9ecef;
        }
        [data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2, [data-testid="stSidebar"] h3 {
            color: var(--text-dark);
        }

        /* --- KNOPPEN STYLING (CORAL) --- */
        
        /* Disabled knop styling (grijs en niet klikbaar) */
        button:disabled {
            background-color: #cccccc !important;
            color: #666666 !important;
            border: 1px solid #999999 !important;
            cursor: not-allowed !important;
        }

        div.stButton > button[kind="primary"]:not(:disabled),
        div.stDownloadButton > button[kind="primary"]:not(:disabled) {
            background-color: var(--primary-color) !important;
            color: white !important;
            border: none;
            font-weight: bold;
            padding: 0.5rem 1rem;
            border-radius: 8px;
            transition: background-color 0.3s;
        }
        
        div.stButton > button[kind="primary"]:hover:not(:disabled),
        div.stDownloadButton > button[kind="primary"]:hover:not(:disabled) {
            background-color: var(--primary-hover) !important;
            color: white !important;
        }

        div.stButton > button, 
        div.stDownloadButton > button {
             border-radius: 8px;
             border: 1px solid var(--primary-color);
             color: var(--text-dark);
        }

        h1.main-title {
            color: var(--text-dark);
            font-weight: 800;
            font-size: 2.5rem;
            margin-bottom: 0;
        }
        .subtitle {
            color: gray;
            font-size: 1.1rem;
            margin-bottom: 2rem;
        }
    </style>
""", unsafe_allow_html=True)

# --- HULPFUNCTIES ---

def get_dxf_center_wgs84(dxf_path):
    """Leest de header van de DXF voor snelle locatiebepaling."""
    try:
        doc = ezdxf.readfile(dxf_path)
        min_x = doc.header.get("$EXTMIN", (0,0,0))[0]
        max_x = doc.header.get("$EXTMAX", (0,0,0))[0]
        min_y = doc.header.get("$EXTMIN", (0,0,0))[1]
        max_y = doc.header.get("$EXTMAX", (0,0,0))[1]
        
        if min_x == 0 and max_x == 0: return None 

        center_x = (min_x + max_x) / 2
        center_y = (min_y + max_y) / 2
        
        transformer = Transformer.from_crs("EPSG:28992", "EPSG:4326", always_xy=True)
        lon, lat = transformer.transform(center_x, center_y)
        return [lat, lon]
    except Exception:
        return None

def generate_map(gpkg_path=None, center_override=None):
    """Genereert de Folium kaart met OSM en data."""
    start_loc = [52.1, 5.1] 
    start_zoom = 8
    
    if center_override:
        start_loc = center_override
        start_zoom = 18 
    elif gpkg_path:
        try:
            layers = fiona.listlayers(gpkg_path)
            if layers:
                for layer in layers:
                    gdf_bounds = gpd.read_file(gpkg_path, layer=layer)
                    if not gdf_bounds.empty:
                        if gdf_bounds.crs is None: gdf_bounds.set_crs("EPSG:28992", inplace=True)
                        gdf_bounds = gdf_bounds.to_crs("EPSG:4326")
                        centroid = gdf_bounds.geometry.unary_union.centroid
                        start_loc = [centroid.y, centroid.x]
                        start_zoom = 18
                        break 
        except: pass 

    m = folium.Map(location=start_loc, zoom_start=start_zoom, tiles="OpenStreetMap")

    if gpkg_path:
        try:
            layers = fiona.listlayers(gpkg_path)
            
            # 20 verschillende kleuren voor maximale variatie
            hex_colors = [
                "#1d2b36", "#e74c3c", "#f39c12",
                "#2c3e50", "#d35400", "#7f8c8d",
                "#27ae60", "#8e44ad", "#2980b9", "#16a085"
            ]
            
            for i, layer_name in enumerate(layers):
                try:
                    gdf = gpd.read_file(gpkg_path, layer=layer_name)
                    if gdf.empty: continue
                    if gdf.crs is None: gdf.set_crs("EPSG:28992", inplace=True)
                    
                    gdf["Laagnaam"] = layer_name
                    
                    if gdf.geom_type.iloc[0] in ['Polygon', 'MultiPolygon']:
                        gdf["Oppervlakte"] = gdf.geometry.area.round(1).astype(str) + " m²"
                    else:
                        gdf["Oppervlakte"] = "-" 

                    gdf_display = gdf.to_crs("EPSG:4326")

                    popup_fields = ["Laagnaam", "Oppervlakte"]
                    popup_aliases = ["Type:", "Oppervlakte:"]
                    
                    for col in gdf_display.columns:
                        if any(x in col.lower() for x in ['materiaal', 'fysiek_voorkomen', 'type']) and col != "Laagnaam":
                            popup_fields.append(col)
                            popup_aliases.append("Info:")
                            break
                    for col in gdf_display.columns:
                         if col.lower() in ['identificatie', 'id']:
                            popup_fields.append(col)
                            popup_aliases.append("ID:")
                            break
                            
                    # Kleur vastzetten met lambda c=...
                    layer_color = hex_colors[i % len(hex_colors)]
                    
                    folium.GeoJson(
                        gdf_display,
                        name=layer_name,
                        style_function=lambda x, c=layer_color: {
                            "fillColor": c,
                            "color": "white", 
                            "weight": 1.5,
                            "fillOpacity": 0.7
                        },
                        marker=folium.CircleMarker(
                            radius=6, 
                            fill_color=layer_color, 
                            color="white", 
                            weight=2, 
                            fill_opacity=1
                        ),
                        tooltip=folium.GeoJsonTooltip(
                            fields=popup_fields[:2], 
                            aliases=popup_aliases[:2],
                            localize=True
                        ),
                        popup=folium.GeoJsonPopup(
                            fields=popup_fields, 
                            aliases=popup_aliases, 
                            localize=True
                        )
                    ).add_to(m)
                    
                except Exception: pass
            folium.LayerControl(position='topright').add_to(m)
        except Exception as e:
            st.error(f"Fout bij laden kaartdata: {e}")

    return m

# --- HOOFD LAYOUT ---

logo_path = "kepa_logo.png"
if os.path.exists(logo_path):
    st.sidebar.image(logo_path, width=220)
else:
    st.sidebar.title("KEPA")

st.sidebar.markdown("---")
st.sidebar.header(" Project Input")

st.markdown('<h1 class="main-title">Revisie Portal</h1>', unsafe_allow_html=True)
st.markdown('<p class="subtitle">Zet meetdata om naar OTL-GIS en CAD tekening.</p>', unsafe_allow_html=True)

# 1. FILE UPLOADER
uploaded_file = st.sidebar.file_uploader("Upload DXF bestand", type=["dxf"])

# 2. METADATA INPUT MET VALIDATIE
form_compleet = False
in_beheerder = ""
in_beheerder_detail = ""
in_opleverdatum = date.today()

if uploaded_file:
    st.sidebar.subheader(" Projectgegevens")
    
    # We voegen een lege string [""] toe aan jouw lijsten voor de placeholder
    beheerders = [""] + ["Gemeente", "Provincie", "Rijkswaterstaat", "Waterschap", "Prorail", "Particulier", "Onbekend", "Overig"]
    beheerders_detail = [""] + ["Amsterdamse Bos", "Havenschap", "Stadsdeel Centrum", "Stadsdeel Nieuw-West", "Stadsdeel Noord", "Stadsdeel Oost", "Stadsdeel West", "Stadsdeel Zuid", "Stadsdeel Zuidoost", "Westpoort", "Weesp", "VVE", "Overig"]
    
    # Selectboxes met format_func om de lege optie netjes te tonen
    in_beheerder = st.sidebar.selectbox(
        "Beheerder *", 
        beheerders,
        format_func=lambda x: "Maak een keuze..." if x == "" else x
    )
    
    in_beheerder_detail = st.sidebar.selectbox(
        "Beheerder (Detail) *", 
        beheerders_detail,
        format_func=lambda x: "Maak een keuze..." if x == "" else x
    )
    
    in_opleverdatum = st.sidebar.date_input("Opleverdatum *", value=date.today())
    
    # Validatie: Check of velden niet leeg zijn
    if in_beheerder != "" and in_beheerder_detail != "":
        form_compleet = True
        # Bereken garantie voor weergave
        garantie_datum = in_opleverdatum.replace(month=in_opleverdatum.month+6) if in_opleverdatum.month <= 6 else in_opleverdatum.replace(year=in_opleverdatum.year+1, month=in_opleverdatum.month-6)
        st.sidebar.info(f"Garantie t/m: {garantie_datum}")
    else:
        st.sidebar.warning("⚠️ Vul de verplichte velden in om te starten.")

if 'map_center' not in st.session_state: st.session_state['map_center'] = None

if uploaded_file:
    if 'last_uploaded' not in st.session_state or st.session_state['last_uploaded'] != uploaded_file.name:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".dxf") as tmp:
            tmp.write(uploaded_file.getbuffer())
            tmp_path = tmp.name
        
        center = get_dxf_center_wgs84(tmp_path)
        if center:
            st.session_state['map_center'] = center
            st.toast("📍 Projectlocatie gevonden!", icon="✅")
        
        st.session_state['last_uploaded'] = uploaded_file.name
        os.remove(tmp_path)

    # KNOP IS DISABLED ALS FORM NIET COMPLEET IS
    if st.sidebar.button(" Start Verwerking", type="primary", use_container_width=True, disabled=not form_compleet):
        
        work_dir = "temp_workspace"
        if os.path.exists(work_dir): shutil.rmtree(work_dir)
        os.makedirs(work_dir)
        
        input_path = os.path.join(work_dir, uploaded_file.name)
        with open(input_path, "wb") as f:
            f.write(uploaded_file.getbuffer())
            
        TEMPLATE_SOURCE = "OTL-Verhardingen_Beheer_2.4.1.gpkg"
        if not os.path.exists(TEMPLATE_SOURCE):
            st.sidebar.error("⚠️ Template GPKG mist op server!")
            st.stop()
        shutil.copy(TEMPLATE_SOURCE, os.path.join(work_dir, TEMPLATE_SOURCE))

        with st.spinner('Bezig met omzetten van CAD naar GIS...'):
            output_dxf_name = uploaded_file.name.replace(".dxf", "_verwerkt.dxf")
            output_dxf_path = os.path.join(work_dir, output_dxf_name)
            
            # Metadata dictionary maken voor de backend
            user_meta = {
                "beheerder": in_beheerder,
                "beheerder_detail": in_beheerder_detail,
                "opleverdatum": in_opleverdatum
            }
            
            try:
                # Metadata doorgeven aan backend
                process_dxf(
                    input_path=input_path,
                    output_path=output_dxf_path,
                    tolerance=0.1,
                    template_path=os.path.join(work_dir, TEMPLATE_SOURCE),
                    user_metadata=user_meta 
                )
                st.sidebar.success("✅ Verwerking klaar!")
            except Exception as e:
                st.sidebar.error(f"Fout tijdens verwerking: {e}")
                st.stop()
        
        st.session_state['gpkg_path'] = None
        st.session_state['dxf_path'] = None
        for f in os.listdir(work_dir):
            if f.endswith(".gpkg") and f != TEMPLATE_SOURCE:
                st.session_state['gpkg_path'] = os.path.join(work_dir, f)
            if f.endswith("_verwerkt.dxf"):
                st.session_state['dxf_path'] = os.path.join(work_dir, f)

if st.session_state.get('gpkg_path') and os.path.exists(st.session_state['gpkg_path']):
    st.sidebar.markdown("---")
    st.sidebar.header(" Resultaten")
    
    with open(st.session_state['gpkg_path'], "rb") as f:
        st.sidebar.download_button(
            label=" OTL-GIS (.gpkg)",
            data=f,
            file_name=os.path.basename(st.session_state['gpkg_path']),
            mime="application/octet-stream",
            type="primary", 
            use_container_width=True
        )
        
    if st.session_state.get('dxf_path') and os.path.exists(st.session_state['dxf_path']):
        with open(st.session_state['dxf_path'], "rb") as f:
            st.sidebar.download_button(
                label=" Hatch Tekening (.dxf)",
                data=f,
                file_name=os.path.basename(st.session_state['dxf_path']),
                mime="application/dxf",
                use_container_width=True
            )

# --- HOOFDSCHERM KAART ---
gpkg_to_show = st.session_state.get('gpkg_path') if st.session_state.get('gpkg_path') and os.path.exists(st.session_state['gpkg_path']) else None
center_to_show = st.session_state.get('map_center')

with st.container(border=True):
    map_obj = generate_map(gpkg_path=gpkg_to_show, center_override=center_to_show)
    st_folium(map_obj, height=750, use_container_width=True)