import streamlit as st
import geopandas as gpd
import pandas as pd
import numpy as np
import os
import zipfile
import tempfile
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from shapely.ops import linemerge, substring
from shapely.geometry import LineString, Polygon
import rasterio
from rasterstats import zonal_stats
import io
import folium
from streamlit_folium import st_folium

# Import untuk ReportLab (PDF)
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import pagesizes
from reportlab.lib.units import inch
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.utils import ImageReader 

# =========================================
# 1. KONFIGURASI HALAMAN UTAMA
# =========================================
st.set_page_config(page_title="ASPAL System", page_icon="🛣️", layout="wide")

# =========================================
# 2. INISIALISASI SESSION STATE GABUNGAN
# =========================================
# State untuk PCI
if 'pci_selesai' not in st.session_state: st.session_state.pci_selesai = False
if 'df_pci' not in st.session_state: st.session_state.df_pci = None
if 'df_detail_pci' not in st.session_state: st.session_state.df_detail_pci = None
if 'peta_bytes_pci' not in st.session_state: st.session_state.peta_bytes_pci = None
if 'grafik_bytes_pci' not in st.session_state: st.session_state.grafik_bytes_pci = None
if 'pdf_bytes_pci' not in st.session_state: st.session_state.pdf_bytes_pci = None
if 'gpkg_bytes_pci' not in st.session_state: st.session_state.gpkg_bytes_pci = None   
if 'seg_gdf_pci' not in st.session_state: st.session_state.seg_gdf_pci = None
if 'excel_bytes_pci' not in st.session_state: st.session_state.excel_bytes_pci = None
if 'master_distress_pci' not in st.session_state: st.session_state.master_distress_pci = None
    
# State untuk SDI
if 'sdi_selesai' not in st.session_state: st.session_state.sdi_selesai = False
if 'df_sdi' not in st.session_state: st.session_state.df_sdi = None
if 'peta_bytes_sdi' not in st.session_state: st.session_state.peta_bytes_sdi = None
if 'grafik_bytes_sdi' not in st.session_state: st.session_state.grafik_bytes_sdi = None
if 'pdf_bytes_sdi' not in st.session_state: st.session_state.pdf_bytes_sdi = None
if 'gpkg_bytes_sdi' not in st.session_state: st.session_state.gpkg_bytes_sdi = None   
if 'seg_gdf_sdi' not in st.session_state: st.session_state.seg_gdf_sdi = None
if 'excel_bytes_sdi' not in st.session_state: st.session_state.excel_bytes_sdi = None
if 'gdf_retak_sdi' not in st.session_state: st.session_state.gdf_retak_sdi = None
if 'gdf_pothole_sdi' not in st.session_state: st.session_state.gdf_pothole_sdi = None
if 'gdf_rutting_sdi' not in st.session_state: st.session_state.gdf_rutting_sdi = None

# ==========================================================
# 3. DATABASE KURVA PCI (STANDAR POLINOMIAL ASTM)
# ==========================================================
DISTRESS_COEFFICIENTS = {
    "alligator_cracking": {"valid_min": 0.1, "valid_max": 100.0, "chart_type": "log", "coefficients": {"low": [11.81030706543641, 14.716555458137659, 5.254969146645571], "medium": [21.641468980402742, 19.850106754347717, 4.129291159183641], "high": [30.698348853111792, 26.819142548270502, 5.653897800825902, -2.0562458038186975]}},
    "bleeding": {"valid_min": 0.1, "valid_max": 100.0, "chart_type": "log", "coefficients": {"low": [0.322105531453021, -0.174525478036764, 1.504981533364469, 1.7947851128512355], "medium": [3.3241213258799234, 4.4914393599717854, 3.3913394399693146, 1.7791635264801178], "high": [5.739963736596728, 7.319479502282341, 7.086578195546586, 3.075398183100246]}},
    "block_cracking": {"valid_min": 0.1, "valid_max": 100.0, "chart_type": "log", "coefficients": {"low": [-0.24506618985589634, 3.5865225127793514, 4.864809264954378], "medium": [2.0044867717721253, 7.656841349922775, 6.197308536458376], "high": [6.054836363905478, 14.23908143225682, 9.26078972426816]}},
    "bumps_and_sags": {"valid_min": 0.1, "valid_max": 10.0, "chart_type": "log", "coefficients": {"low": [7.432162927047138, 13.039704473951343, 12.681378093050206, 5.724727901683913], "medium": [24.0510616920146, 24.990801854094467, 17.775330890078287, 10.801916366371941], "high": [53.60256635094115, 38.448582442650476, 6.043250542072165]}},
    "corrugation": {"valid_min": 0.1, "valid_max": 100.0, "chart_type": "log", "coefficients": {"low": [2.079645270606063, 6.171432244797827, 6.315591279321682], "medium": [16.00629930412681, 17.42721794385406, 5.877793175640937], "high": [33.3926340071105, 25.16870941731629, 2.8124490453717783]}},
    "depression": {"valid_min": 0.1, "valid_max": 100.0, "chart_type": "log", "coefficients": {"low": [3.289661412212115, 0.6252836980242593, 10.964451357494776, 6.985101871322738, -3.5060456381316634], "medium": [7.614400618717246, 3.769091495110122, 15.559272647446328, 7.4759184063961435, -4.844857281255507], "high": [15.948684201140967, 9.394030228694518, 15.510038701741792, 6.024061067798476, -4.422031207577385]}},
    "edge_cracking": {"valid_min": 0.1, "valid_max": 20.0, "chart_type": "log", "coefficients": {"low": [3.211409886900335, 4.962171357579026, 2.4025247047022504], "medium": [9.38292671744566, 9.608124129013511, 4.3396040463768015], "high": [15.664668848033116, 15.377839697997818, 7.293639747414968]}},
    "joint_reflection_cracking": {"valid_min": 0.1, "valid_max": 30.0, "chart_type": "log", "coefficients": {"low": [2.599637864741199, 7.548320214895041, 6.035620574755107], "medium": [7.381829148226155, 13.409513965059839, 14.487037415701224, 2.604980254156252, -4.850878879939305], "high": [15.708313735606989, 21.751135896400438, 22.35224402708013, 12.321810537440612, -6.404620335343928, -4.635950974262599]}},
    "lane_shoulder_drop_off": {"valid_min": 0.5, "valid_max": 15.0, "chart_type": "log", "coefficients": {"low": [1.7330434661503789, 2.2479876898054147, 8.492414136167602], "medium": [3.1018055744999637, 0.9748217115112947, 15.880585262031168], "high": [5.588313880268777, 5.140183474398842, 22.853461108687]}},
    "longitudinal_transverse_cracking": {"valid_min": 0.1, "valid_max": 30.0, "chart_type": "log", "coefficients": {"low": [2.1831950830355784, 8.216281290902756, 6.454535216997481], "medium": [8.67907975505781, 15.31510859022643, 6.266899017617167], "high": [17.723256518500968, 24.494224894065376, 19.119102009732657, 4.182799471044728, -4.54528270804704]}},
    "patching_and_utility_cut_patching": {"valid_min": 0.1, "valid_max": 50.0, "chart_type": "log", "coefficients": {"low": [1.8384973325128793, 7.632382033226324, 6.507032451977727], "medium": [8.929258996532853, 14.056348057074333, 8.632199015982204], "high": [18.05379050800971, 18.570848940912526, 15.195881745184764, 4.634730303593308, -4.282412070934171]}},
    "polished_aggregate": {"valid_min": 0.1, "valid_max": 100.0, "chart_type": "log", "coefficients": {"low": [0.16789523734739448, -0.1540196727901293, 1.3119789947979732, 1.8893029190226103], "medium": [0.16789523734739448, -0.1540196727901293, 1.3119789947979732, 1.8893029190226103], "high": [0.16789523734739448, -0.1540196727901293, 1.3119789947979732, 1.8893029190226103]}},
    "potholes": {"valid_min": 0.01, "valid_max": 10.0, "chart_type": "log", "coefficients": {"low": [58.574736886209024, 41.33443795417745, 2.307796828048822, -2.1009955078236295], "medium": [91.64010353001909, 65.40209466565399, 5.262639941411422, -3.0315939096971647], "high": [109.3323221736685, 56.29275767507379, -0.3934988005267144, -3.0772594822340906]}},
    "railroad_crossing": {"valid_min": 1.0, "valid_max": 40.0, "chart_type": "log", "coefficients": {"low": [0.5029257285845483, 5.749904110236876, 3.9821246564884927], "medium": [5.254000147362391, 11.422544085878982, 37.11519773645627, -16.910805393020667], "high": [19.340445737998966, 18.23862361469355, 53.47893872454435, -25.833753655092053]}},
    "raveling": {"valid_min": 0.1, "valid_max": 100.0, "chart_type": "log", "coefficients": {"low": [1.3364599508023929, 2.5381741896122354, 2.2521180739800015], "medium": [8.547489838241672, 5.50999488268751, 2.8976928296955693, 1.6574499250743635], "high": [15.287212880590612, 12.713212203226831, 11.470876324575428, 5.038585870016526, -3.053732169374685]}},
    "rutting": {"valid_min": 0.1, "valid_max": 100.0, "chart_type": "log", "coefficients": {"low": [7.166892989037484, 15.480212329682374, 7.204541988003889, -2.056356612335577], "medium": [17.07939650585985, 23.042120387468632, 7.267187817332259, -2.992255156045317], "high": [26.55168225206831, 25.25096737467123, 9.616695504985966, 2.709242740203141, -2.9515698599184823]}},
    "shoving": {"valid_min": 0.1, "valid_max": 50.0, "chart_type": "log", "coefficients": {"low": [4.413031420403843, 10.016571238629226, 5.535372489599174], "medium": [9.294794348392088, 16.13985140155734, 10.13135831372187], "high": [17.76171628163815, 19.396982653439373, 15.785981300450224, 3.016185582467493, -3.721051791249142]}},
    "slippage_cracking": {"valid_min": 0.1, "valid_max": 100.0, "chart_type": "log", "coefficients": {"low": [4.248092940884805, 14.667313248550373, 9.293396804566436, -2.1419610585945215], "medium": [10.838179259718574, 20.247013730206078, 13.130140952344512, -0.24505223183555103, -2.039447457156605], "high": [17.75416989366112, 32.470845743485654, 29.87775942697452, -5.259819524132407, -13.039225691175872, 4.346079607761759]}},
    "swell": {"valid_min": 1.0, "valid_max": 30.0, "chart_type": "log", "coefficients": {"low": [0.2776553410390519, 12.183535319307815], "medium": [9.7262404794108, 26.04919313918952], "high": [33.021164499056795, 9.27090451459351, 10.656396250845475]}},
    "weathering": {"valid_min": 0.1, "valid_max": 100.0, "chart_type": "log", "coefficients": {"low": [1.3364599508023929, 2.5381741896122354, 2.2521180739800015], "medium": [8.547489838241672, 5.50999488268751, 2.8976928296955693, 1.6574499250743635], "high": [15.287212880590612, 12.713212203226831, 11.470876324575428, 5.038585870016526, -3.053732169374685]}}
}

CDV_ASPHALT_COEFFICIENTS = {
    "q1": [0.0, 0.9999999999999999],
    "q2": [-4.767827657379179, 0.9204018317853505, -0.001751870485036131],
    "q3": [-7.225000000000371, 0.8352547729618199, -0.0013491357069143565],
    "q4": [-12.07301341589303, 0.8359713622291058, -0.0013966718266254004],
    "q5": [-12.825902992776388, 0.7672961816305504, -0.0011578947368421164],
    "q6": [-15.061893704850672, 0.7418633900928823, -0.0010566950464396384],
    "q7": [-18.186068111455448, 0.8352115583075376, -0.0016336429308565648],
    "q8": [-11.661042311659124, 0.5804493464050489, 0.0012355521155837934, -9.988820089439906e-06],
    "q9": [-10.821138630889415, 0.5307870370368919, 0.0020743464052291917, -1.4092420593968467e-05],
    "q10": [-33.76740887897169, 2.045383029591749, -0.033305999495864574, 0.00035829181520672087, -1.7915566098874656e-06, 3.173417753551723e-09]
}

# ==========================================================
# 4. FUNGSI SPASIAL & MATEMATIKA (SHARED)
# ==========================================================
def read_zip_shapefile(uploaded_file, tmpdir):
    zip_path = os.path.join(tmpdir, uploaded_file.name)
    with open(zip_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    extract_dir = os.path.join(tmpdir, uploaded_file.name.replace('.zip', ''))
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_dir)
    for root, dirs, files in os.walk(extract_dir):
        for file in files:
            if file.endswith(".shp"):
                return gpd.read_file(os.path.join(root, file))
    return None

def hitung_width(gdf):
    gdf = gdf.copy()
    width_list = []
    for geom in gdf.geometry:
        if geom.is_empty:
            width_list.append(np.nan)
            continue
        mrr = geom.minimum_rotated_rectangle
        coords = list(mrr.exterior.coords)
        edges = [LineString([coords[i], coords[i+1]]).length for i in range(4)]
        width_list.append(min(edges) * 1000)
    gdf["WIDTH_MM"] = width_list
    return gdf

def hitung_diameter_pothole(gdf):
    gdf = gdf.copy()
    diameter_list = []
    for geom in gdf.geometry:
        if geom.is_empty:
            diameter_list.append(np.nan)
            continue
        mrr = geom.minimum_rotated_rectangle
        coords = list(mrr.exterior.coords)
        edges = [LineString([coords[i], coords[i+1]]).length for i in range(4)]
        unique_edges = sorted(list(set([round(e,5) for e in edges])))
        if len(unique_edges) >= 2:
            diameter = (max(unique_edges) + min(unique_edges)) / 2
        else:
            diameter = unique_edges[0]
        diameter_list.append(diameter * 1000)
    gdf["DIAMETER_MM"] = diameter_list
    return gdf

def hitung_depth(gdf, dsm_path, buffer_distance=0.3):
    """Menghitung kedalaman dalam milimeter (mm) untuk PCI"""
    with rasterio.open(dsm_path) as DSM:
        dsm_crs = DSM.crs
        nodata_val = DSM.nodata
        
    buffer_outer = gdf.geometry.buffer(buffer_distance)
    ring_geom = buffer_outer.difference(gdf.geometry)
    hole_geom_dsm = gdf.geometry.to_crs(dsm_crs)
    ring_geom_dsm = ring_geom.to_crs(dsm_crs)
    
    stats_hole = zonal_stats(hole_geom_dsm, dsm_path, stats=["median"], nodata=nodata_val, all_touched=True)
    stats_ring = zonal_stats(ring_geom_dsm, dsm_path, stats=["median"], nodata=nodata_val, all_touched=True)
    
    depth_list = []
    for i in range(len(gdf)):
        z_hole = stats_hole[i]["median"]
        z_ref = stats_ring[i]["median"]
        
        if z_hole is not None and z_ref is not None:
            depth = (z_ref - z_hole) * 1000
            depth = max(0, min(depth, 80))
        else:
            depth = 0
            
        depth_list.append(depth)
        
    gdf = gdf.copy()
    gdf["DEPTH_MM"] = depth_list
    return gdf


def hitung_depth_cm(gdf, dsm_path, buffer_distance=0.3):
    """Menghitung kedalaman dalam centimeter (cm) untuk SDI"""
    with rasterio.open(dsm_path) as DSM:
        dsm_crs = DSM.crs
        nodata_val = DSM.nodata
        
    buffer_outer = gdf.geometry.buffer(buffer_distance)
    ring_geom = buffer_outer.difference(gdf.geometry)
    hole_geom_dsm = gdf.geometry.to_crs(dsm_crs)
    ring_geom_dsm = ring_geom.to_crs(dsm_crs)
    
    stats_hole = zonal_stats(hole_geom_dsm, dsm_path, stats=["median"], nodata=nodata_val, all_touched=True)
    stats_ring = zonal_stats(ring_geom_dsm, dsm_path, stats=["median"], nodata=nodata_val, all_touched=True)
    
    depth_list = []
    for i in range(len(gdf)):
        z_hole = stats_hole[i]["median"]
        z_ref = stats_ring[i]["median"]
        
        if z_hole is not None and z_ref is not None:
            depth = (z_ref - z_hole) * 100
            depth = max(0, min(depth, 15))
        else:
            depth = 0
            
        depth_list.append(depth)
        
    gdf = gdf.copy()
    gdf["kedalaman_calc"] = depth_list
    return gdf

def tentukan_severity(distress_type, row):
    distress = distress_type.lower()
    depth = row.get("Depth_mm", 0)
    width = row.get("Width_mm", 0)
    diameter = row.get("Diameter_mm", 0)
    area = row.geometry.area if row.geometry is not None else 0

    if "alligator" in distress: return "Low" if width < 10 else "Medium" if width <= 25 else "High"
    if "bleeding" in distress: return "Low" if area < 2 else "Medium" if area < 6 else "High"
    if "block" in distress: return "Low" if width < 10 else "Medium" if width <= 25 else "High"
    if "bump" in distress or "sag" in distress: return "Low" if depth <= 10 else "Medium" if depth <= 25 else "High"
    if "corrugation" in distress: return "Low" if depth <= 10 else "Medium" if depth <= 25 else "High"
    if "depression" in distress: return "Low" if depth <= 25 else "Medium" if depth <= 50 else "High"
    if "edge" in distress: return "Low" if width < 10 else "Medium" if width <= 25 else "High"
    if "joint" in distress or "reflection" in distress: return "Low" if width < 10 else "Medium" if width <= 75 else "High"
    if "shoulder" in distress: return "Low" if depth <= 50 else "Medium" if depth <= 100 else "High"
    if "longitudinal" in distress or "transverse" in distress: return "Low" if width < 10 else "Medium" if width <= 75 else "High"
    if "patch" in distress: return "Low" if area < 1 else "Medium" if area < 3 else "High"
    if "polished" in distress: return "Low" if area < 5 else "Medium" if area < 15 else "High"
    if "pothole" in distress:
        if depth <= 25: return "Low" if diameter < 450 else "Medium"
        elif depth <= 50: return "Low" if diameter < 200 else "Medium" if diameter < 450 else "High"
        else: return "Medium" if diameter < 450 else "High"
    if "railroad" in distress: return "Low" if depth <= 25 else "Medium" if depth <= 75 else "High"
    if "ravel" in distress: return "Low" if area < 5 else "Medium" if area < 20 else "High"
    if "rutting" in distress: return "Low" if depth <= 13 else "Medium" if depth <= 25 else "High"
    if "shoving" in distress: return "Low" if depth <= 10 else "Medium" if depth <= 25 else "High"
    if "slippage" in distress: return "Low" if width < 10 else "Medium" if width <= 40 else "High"
    if "swell" in distress: return "Low" if depth <= 25 else "Medium" if depth <= 75 else "High"
    if "weather" in distress: return "Low" if area < 5 else "Medium" if area < 20 else "High"
    return "Low"

def evaluate_polynomial(coeffs, x):
    return sum(c * (x ** i) for i, c in enumerate(coeffs))

def lookup_dv(distress_type, severity, density):
    key = distress_type.lower().replace(" ", "_").replace("-", "_")
    if key not in DISTRESS_COEFFICIENTS: return 0.0
    data = DISTRESS_COEFFICIENTS[key]
    coeffs = data["coefficients"].get(severity.lower())
    if not coeffs: return 0.0
    density = max(data["valid_min"], min(data["valid_max"], density))
    val_to_eval = math.log10(density) if data["chart_type"] == "log" else density
    dv = evaluate_polynomial(coeffs, val_to_eval)
    return max(0.0, min(100.0, dv))

def lookup_cdv_asphalt(q, total_deduct_value):
    q_lookup = f"q{min(max(int(q), 1), 10)}"
    coeffs = CDV_ASPHALT_COEFFICIENTS.get(q_lookup)
    tdv_clamped = max(0.0, min(200.0, total_deduct_value))
    cdv = evaluate_polynomial(coeffs, tdv_clamped)
    return max(0.0, min(100.0, cdv))

def rating_pci(pci):
    if pci > 85: return "Good"
    elif pci > 70: return "Satisfactory"
    elif pci > 55: return "Fair"
    elif pci > 40: return "Poor"
    elif pci > 25: return "Very Poor"
    elif pci > 10: return "Serious"
    else: return "Failed"

def hitung_sdi(persen_retak, lebar_retak, jumlah_lubang, kedalaman_rutting):
    if persen_retak == 0: sdi1 = 0
    elif persen_retak < 10: sdi1 = 5
    elif persen_retak <= 30: sdi1 = 20
    else: sdi1 = 40

    sdi2 = sdi1 * 2 if lebar_retak > 3 else sdi1

    if jumlah_lubang == 0: sdi3 = sdi2
    elif jumlah_lubang < 10: sdi3 = sdi2 + 15
    elif jumlah_lubang <= 50: sdi3 = sdi2 + 75
    else: sdi3 = sdi2 + 225

    if kedalaman_rutting == 0: sdi4 = sdi3
    elif kedalaman_rutting < 1: sdi4 = sdi3 + (5 * 0.5)
    elif kedalaman_rutting <= 3: sdi4 = sdi3 + (5 * 2)
    else: sdi4 = sdi3 + (5 * 4)

    if sdi4 < 50: kondisi = "Baik"
    elif sdi4 <= 100: kondisi = "Sedang"
    elif sdi4 <= 150: kondisi = "Rusak Ringan"
    else: kondisi = "Rusak Berat"
        
    return sdi1, sdi2, sdi3, sdi4, kondisi

def metric_card(label, value, value_color="#4da6ff", bg_color="#1E2A38", text_color="#cbd5e1"):
    return f"""
    <div style="background-color: {bg_color}; padding: 15px; border-radius: 8px; border: 1px solid #2d3e50; text-align: center; height: 100%;">
        <p style="margin: 0px; font-size: 14px; color: {text_color};">{label}</p>
        <h2 style="margin: 5px 0px 0px 0px; color: {value_color}; font-size: 24px; font-weight: bold;">{value}</h2>
    </div>
    """

# =========================================
# 5. SIDEBAR & NAVIGASI SISTEM ASPAL
# =========================================
with st.sidebar:
    st.title("🛣️ ASPAL")
    st.caption("Analisis Spasial Perkerasan Jalan")
    
    menu = st.radio("Pilih Modul Analisis:", [
        "🏠 Beranda", 
        "📈 Modul PCI (Pavement Condition Index)", 
        "📉 Modul SDI (Surface Distress Index)",
        "📊 Komparasi (PCI vs SDI)"
    ])
    
    st.divider()
    st.header("📝 Informasi Umum Survey")
    lokasi = st.text_input("Lokasi Survey", "Jl. Contoh Raya")
    sta_umum = st.text_input("STA Umum", "0+000 - 1+000")
    surveyor = st.text_input("Nama Surveyor", "Nama Anda")
    tanggal = st.text_input("Tanggal Survey", "28 Februari 2026")
    instansi = st.text_input("Instansi", "Universitas Diponegoro")
    
    st.header("⚙️ Parameter Jalan")
    lebar_jalan = st.number_input("Lebar Jalan (m)", value=3.5, step=0.1)
    interval_segmen = st.number_input("Interval Segmen (m)", value=100, step=10)
    epsg_code = st.number_input("Kode EPSG UTM Lokal (Misal: 32749)", value=32749, step=1)
   
    st.divider()
    if st.button("🔄 Reset Sistem ASPAL", use_container_width=True):
        st.session_state.clear()
        st.rerun()

# =========================================
# 6. ROUTING HALAMAN BERDASARKAN MENU
# =========================================

if menu == "🏠 Beranda":
    # Hero Section
    st.markdown("<h1 style='text-align: center; color: #4da6ff;'>🛣️ Sistem ASPAL</h1>", unsafe_allow_html=True)
    st.markdown("<h4 style='text-align: center; color: #cbd5e1;'>Analisis Spasial Perkerasan Jalan Berbasis GIS</h4>", unsafe_allow_html=True)
    st.divider()

    st.markdown("""
    Selamat datang di **ASPAL**, platform cerdas untuk otomatisasi evaluasi kondisi infrastruktur jalan. 
    Sistem ini mengintegrasikan dua metode penilaian standar (PCI dan SDI) dengan analisis geospasial untuk menghasilkan rekomendasi pemeliharaan yang akurat, cepat, dan terukur.
    """)

    st.markdown("<br>", unsafe_allow_html=True)

    # ==========================================
    # KARTU ALUR KERJA (WORKFLOW)
    # ==========================================
    st.markdown("### 🚀 Alur Kerja Sistem")
    col_w1, col_w2, col_w3 = st.columns(3)
    
    with col_w1:
        st.markdown("""
        <div style='background-color: #1E2A38; padding: 25px 20px; border-radius: 10px; border: 1px solid #2d3e50; text-align: center; min-height: 210px;'>
            <h1 style='margin: 0 0 10px 0; font-size: 40px;'>📂</h1>
            <h4 style='color: white; margin: 0 0 15px 0;'>1. Input Data</h4>
            <p style='font-size: 13.5px; color: #cbd5e1; line-height: 1.5; margin: 0;'>Unggah Shapefile Jalan, layer polygon kerusakan (Retak, Lubang, dll), dan data elevasi DSM.</p>
        </div>
        """, unsafe_allow_html=True)
        
    with col_w2:
        st.markdown("""
        <div style='background-color: #1E2A38; padding: 25px 20px; border-radius: 10px; border: 1px solid #2d3e50; text-align: center; min-height: 210px;'>
            <h1 style='margin: 0 0 10px 0; font-size: 40px;'>⚙️</h1>
            <h4 style='color: white; margin: 0 0 15px 0;'>2. Proses Spasial</h4>
            <p style='font-size: 13.5px; color: #cbd5e1; line-height: 1.5; margin: 0;'>Sistem otomatis memotong segmen, menghitung luas atau kedalaman kerusakan, dan mengkalkulasi indeks.</p>
        </div>
        """, unsafe_allow_html=True)
        
    with col_w3:
        st.markdown("""
        <div style='background-color: #1E2A38; padding: 25px 20px; border-radius: 10px; border: 1px solid #2d3e50; text-align: center; min-height: 210px;'>
            <h1 style='margin: 0 0 10px 0; font-size: 40px;'>📊</h1>
            <h4 style='color: white; margin: 0 0 15px 0;'>3. Laporan & Peta</h4>
            <p style='font-size: 13.5px; color: #cbd5e1; line-height: 1.5; margin: 0;'>Unduh hasil akhir dalam bentuk PDF berstandar ASTM, Peta Spasial (GPKG), dan tabel Excel.</p>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("<br><br>", unsafe_allow_html=True)

    # ==========================================
    # PENJELASAN PARAMETER
    # ==========================================
    st.markdown("### 🔍 Parameter yang Dievaluasi")
    
    col_b1, col_b2 = st.columns(2)
    
    with col_b1:
        st.info("**Pavement Condition Index (PCI)**")
        st.markdown("""
        Mengkalkulasi tingkat kerusakan kompleks berdasarkan **19 jenis kerusakan** (Standar ASTM D6433):
        
        * **Retak:** Kulit Buaya, Blok, Pinggir, Memanjang/Melintang, Refleksi, Selip.
        * **Deformasi:** Alur (*Rutting*), Keriting, Amblesan, Jembulan, Sungkur.
        * **Permukaan Dasar:** Lubang, Kegemukan Aspal, Pelepasan Butir, Pelapukan.
        * **Lain-lain:** Penurunan Bahu Jalan, Tambalan, Perlintasan KA.
        """)
        
    with col_b2:
        st.success("**Surface Distress Index (SDI)**")
        st.markdown("""
        Metode fungsional dari Bina Marga yang berfokus pada **4 indikator utama** kerusakan visual:
        
        * **Luas Retak:** Persentase area retak terhadap total luas segmen.
        * **Lebar Retak:** Rata-rata bukaan celah retak (< 1mm, 1-3mm, > 3mm).
        * **Jumlah Lubang:** Total titik lubang (*potholes*) per segmen.
        * **Kedalaman Alur:** Penurunan jejak roda (*rutting*) yang diekstrak secara 3D dari DSM.
        """)
        
    st.divider()
    
    # CTA Formal
    st.info("**Panduan Penggunaan:** Untuk memulai tahapan analisis dan evaluasi kondisi perkerasan jalan, silakan memilih **Modul PCI** atau **Modul SDI** melalui panel navigasi yang tersedia di sebelah kiri.")

# =========================================
# MODUL PCI
# =========================================
elif menu == "📈 Modul PCI (Pavement Condition Index)":
    st.title("🛣️ Modul Evaluasi PCI")
    st.markdown("Otomatisasi perhitungan Pavement Condition Index (PCI) metode ASTM D6433.")
    
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("📁 1. Data Dasar")
        jalan_file = st.file_uploader("Upload Shapefile Jalan (.zip)", type="zip", key="pci_jalan")
        dsm_mode = st.radio("Cara Input Data DSM:", ["Upload File .tif", "Paste Link Google Drive"], key="pci_dsm_mode")
        dsm_file = None; dsm_link = ""
        if dsm_mode == "Upload File .tif":
            dsm_file = st.file_uploader("Upload Data DSM (.tif)", type="tif", key="pci_dsm_file")
        else:
            dsm_link = st.text_input("Paste Link Shareable Google Drive (.tif)", key="pci_dsm_link")
            
    with col2:
        st.subheader("⚠️ 2. Data Kerusakan (Distress)")
        distress_keys = list(DISTRESS_COEFFICIENTS.keys())
        distress_options = [k.replace('_', ' ').title() for k in distress_keys]
        selected_distress = st.multiselect("Pilih Kerusakan yang Ditemukan:", distress_options)
        
        uploaded_distress = {}
        for d in selected_distress:
            file = st.file_uploader(f"Upload SHP {d} (.zip)", type="zip", key=f"pci_{d}")
            if file:
                original_key = d.lower().replace(' ', '_')
                uploaded_distress[original_key] = file

    if st.button("🚀 Proses & Hitung PCI", type="primary", use_container_width=True):
        is_dsm_valid = (dsm_mode == "Upload File .tif" and dsm_file is not None) or (dsm_mode == "Paste Link Google Drive" and dsm_link != "")
        if not jalan_file or not uploaded_distress or not is_dsm_valid:
            st.error("⚠️ Mohon lengkapi Shapefile Jalan, Data DSM, dan minimal 1 Data Kerusakan.")
        else:
            with st.spinner("Memproses Analisis Geospasial PCI..."):
                with tempfile.TemporaryDirectory() as tmpdir:
                    try:
                        # 1. BACA JALAN
                        jalan = read_zip_shapefile(jalan_file, tmpdir)
                        if jalan.crs is None: st.error("CRS shapefile jalan tidak terdefinisi!"); st.stop()
                        if jalan.crs.to_epsg() != epsg_code: jalan = jalan.to_crs(epsg=epsg_code)
                        union_geom = jalan.geometry.union_all()
                        merged_line = linemerge(union_geom) if union_geom.geom_type == "MultiLineString" else union_geom
                        panjang_total = merged_line.length
                        segments = [substring(merged_line, start, min(start + interval_segmen, panjang_total)) 
                                    for start in np.arange(0, panjang_total, interval_segmen)]
                        seg_gdf = gpd.GeoDataFrame(geometry=segments, crs=jalan.crs)
                        seg_gdf["Segmen"] = range(1, len(seg_gdf)+1)
                        seg_gdf["geometry"] = seg_gdf.buffer(lebar_jalan / 2, cap_style=2)
                        seg_gdf["Unit_Area"] = seg_gdf.geometry.area
                        
                        # 2. DSM
                        dsm_path = os.path.join(tmpdir, "dsm.tif")
                        if dsm_mode == "Upload File .tif":
                            with open(dsm_path, "wb") as f: f.write(dsm_file.getbuffer())
                        else:
                            import gdown, re
                            match = re.search(r"/d/([a-zA-Z0-9_-]+)", dsm_link)
                            if match: gdown.download(id=match.group(1), output=dsm_path, quiet=False)
                            else: st.error("❌ Link Google Drive tidak valid."); st.stop()

                        # 3. PROSES DISTRESS
                        distress_layers = {}
                        for key, file in uploaded_distress.items():
                            gdf = read_zip_shapefile(file, tmpdir)
                            if gdf is not None: distress_layers[key] = gdf

                        all_distress_list = []
                        target_crs = seg_gdf.crs
                        for nama_distress, gdf in distress_layers.items():
                            if gdf.empty: continue
                            gdf = gdf.copy()
                            if gdf.crs is None: gdf.set_crs(target_crs, inplace=True, allow_override=True)
                            elif gdf.crs != target_crs: gdf = gdf.to_crs(target_crs)
                            
                            if any(x in nama_distress for x in ["crack", "alligator", "block", "long", "slip", "joint", "edge"]):
                                gdf = hitung_width(gdf)
                            if any(x in nama_distress for x in ["rutting", "depression", "corrugation", "bump", "sag", "shoving", "swell", "shoulder", "railroad"]):
                                gdf = hitung_depth(gdf, dsm_path)
                            if "pothole" in nama_distress:
                                gdf = hitung_depth(gdf, dsm_path)
                                gdf = hitung_diameter_pothole(gdf)
                                
                            gdf["Distress_Type"] = nama_distress
                            gdf["Severity"] = gdf.apply(lambda row: tentukan_severity(nama_distress, row), axis=1)
                            gdf["Priority"] = gdf["Severity"].map({"High": 3, "Medium": 2, "Low": 1}).fillna(1)
                            all_distress_list.append(gdf)

                        # 4. OVERLAY & FLATTEN
                        df_detail_list = []
                        if all_distress_list:
                            cleaned_layers = []
                            for gdf in all_distress_list:
                                if gdf.crs != seg_gdf.crs: gdf = gdf.to_crs(seg_gdf.crs)
                                gdf = gdf.set_crs(seg_gdf.crs, allow_override=True)
                                cleaned_layers.append(gdf)
                            master_distress = gpd.GeoDataFrame(pd.concat(cleaned_layers, ignore_index=True), crs=seg_gdf.crs)
                            master_distress = master_distress.sort_values(by="Priority", ascending=False).reset_index(drop=True)
                            st.session_state.master_distress_pci = master_distress.copy()
                        
                            accumulated_geom = Polygon()
                            cleaned_geometries = []
                            for idx, row in master_distress.iterrows():
                                geom = row.geometry
                                new_geom = geom.difference(accumulated_geom) if not geom.is_empty else geom
                                cleaned_geometries.append(new_geom)
                                if not new_geom.is_empty: accumulated_geom = accumulated_geom.union(new_geom)
                            master_distress["geometry"] = cleaned_geometries
                            master_distress = master_distress[~master_distress.geometry.is_empty]
                            inter_all = gpd.overlay(master_distress, seg_gdf, how="intersection")
                        
                            if not inter_all.empty:
                                inter_all["Area_Intersect"] = inter_all.geometry.area
                                agg_df = inter_all.groupby(['Segmen', 'Distress_Type', 'Severity', 'Unit_Area'])['Area_Intersect'].sum().reset_index()
                                for _, row in agg_df.iterrows():
                                    density = max(0, min(100, (row['Area_Intersect'] / row['Unit_Area']) * 100))
                                    dv = lookup_dv(row['Distress_Type'], row['Severity'], density)
                                    df_detail_list.append({"Segmen": row["Segmen"], "Distress": row['Distress_Type'], "Severity": row["Severity"], "Density": density, "DV": dv})

                        df_detail = pd.DataFrame(df_detail_list)
                        
                        # 5. HITUNG PCI
                        df_pci_list = []
                        for seg_id in seg_gdf["Segmen"]:
                            df_seg = df_detail[df_detail["Segmen"] == seg_id] if not df_detail.empty else pd.DataFrame()
                            if df_seg.empty:
                                df_pci_list.append({"Segmen": seg_id, "TDV": 0, "q": 0, "CDV": 0, "PCI": 100})
                                continue
                            
                            dvs = np.array(sorted(df_seg["DV"].tolist(), reverse=True))
                            hdv = dvs[0]
                            m = min(10.0, 1 + (9.0 / 95.0) * (100.0 - hdv))
                            m_int = int(np.floor(m))
                            entered_dvs = dvs[:m_int + 1].copy() if len(dvs) > m_int else np.copy(dvs)
                            if len(dvs) > m_int: entered_dvs[-1] *= (m - m_int)
                            
                            max_cdv = 0.0
                            while True:
                                q = np.sum(entered_dvs > 2.0)
                                if q == 0 and max_cdv != 0: break
                                cdv_current = lookup_cdv_asphalt(q, np.sum(entered_dvs))
                                if cdv_current > max_cdv: max_cdv = cdv_current
                                if q <= 1: break
                                gt_2_indices = np.where(entered_dvs > 2.0)[0]
                                if len(gt_2_indices) > 0: entered_dvs[gt_2_indices[np.argmin(entered_dvs[gt_2_indices])]] = 2.0
                                else: break
                                
                            df_pci_list.append({"Segmen": seg_id, "TDV": round(np.sum(entered_dvs), 2), "q": q, "CDV": round(max_cdv, 2), "PCI": round(max(0, min(100, 100 - max_cdv)), 2)})
                        
                        df_pci = pd.DataFrame(df_pci_list)
                        df_pci["Rating"] = df_pci["PCI"].apply(rating_pci)
                        df_pci["STA"] = df_pci["Segmen"].apply(lambda x: f"{(x-1)*interval_segmen} - {x*interval_segmen} m")
                        
                        seg_gdf = seg_gdf.merge(df_pci, on="Segmen", how="left")
                        seg_gdf["PCI"] = seg_gdf["PCI"].fillna(100)
                        seg_gdf["Rating"] = seg_gdf["Rating"].fillna("Good")

                        # Peta & Visualisasi Matplotlib
                        warna_pci = {"Good": "#006400", "Satisfactory": "#8FBC8F", "Fair": "#FFFF00", "Poor": "#FF6347", "Very Poor": "#FF4500", "Serious": "#8B0000", "Failed": "#A9A9A9"}
                        fig_map, ax_map = plt.subplots(figsize=(10,6))
                        seg_plot = seg_gdf.copy()
                        seg_plot["geometry"] = seg_plot.geometry.buffer(4)
                        legend_handles = []
                        for rating, warna in warna_pci.items():
                            subset = seg_plot[seg_plot["Rating"] == rating]
                            if not subset.empty: 
                                subset.plot(ax=ax_map, color=warna, edgecolor="black", label=f"{rating}")
                                legend_handles.append(mpatches.Patch(color=warna, label=f"{rating} ({len(subset)})"))
                        for idx, row in seg_gdf.iterrows():
                            centroid = row.geometry.centroid
                            ax_map.text(centroid.x, centroid.y, f"S{row['Segmen']}\n{row['PCI']:.0f}", fontsize=7, weight="bold", ha="center", va="center", bbox=dict(facecolor="white", alpha=0.8, boxstyle="round,pad=0.2", edgecolor="gray", lw=0.5))
                        if legend_handles: ax_map.legend(handles=legend_handles, loc="best", title="Kategori PCI", fontsize=8)
                        ax_map.axis("off")
                        peta_path = os.path.join(tmpdir, "peta_pci.png")
                        plt.savefig(peta_path, dpi=300, bbox_inches='tight')
                        plt.close(fig_map)
                        
                        fig_bar, ax_bar = plt.subplots(figsize=(6,4))
                        rekap = seg_gdf["Rating"].value_counts()
                        rekap.plot(kind="bar", color=[warna_pci.get(x, "grey") for x in rekap.index], edgecolor="black", ax=ax_bar)
                        plt.xticks(rotation=45)
                        plt.tight_layout()
                        grafik_path = os.path.join(tmpdir, "grafik_pci.png")
                        plt.savefig(grafik_path, dpi=300)
                        plt.close(fig_bar)
                        
                        # Generate Full PDF Report (PCI)
                        pdf_path = os.path.join(tmpdir, "Laporan_PCI.pdf")
                        doc = SimpleDocTemplate(pdf_path, pagesize=pagesizes.landscape(pagesizes.A4), rightMargin=30, leftMargin=30, topMargin=30, bottomMargin=30)
                        elements = []
                        styles = getSampleStyleSheet()
                        cover_style = ParagraphStyle('cover', parent=styles['Title'], alignment=TA_CENTER)
                        header_style = ParagraphStyle('header', parent=styles['Normal'], alignment=TA_LEFT, fontSize=12, spaceAfter=10, textColor=colors.HexColor("#1f2937"))
                        
                        total_area = seg_gdf["Unit_Area"].sum() if len(seg_gdf) > 0 else 0
                        rata_pci = round((seg_gdf["PCI"] * seg_gdf["Unit_Area"]).sum() / total_area, 2) if total_area > 0 else 0
                        kondisi_dominan = seg_gdf["Rating"].value_counts().idxmax() if not seg_gdf["Rating"].empty else "-"

                        elements.append(Paragraph(instansi, cover_style)); elements.append(Spacer(1, 0.3*inch))
                        elements.append(Paragraph("LAPORAN SURVEY", cover_style)); elements.append(Spacer(1, 0.3*inch))
                        elements.append(Paragraph("PAVEMENT CONDITION INDEX (PCI)", cover_style)); elements.append(Spacer(1, 1*inch))
                        elements.append(Paragraph(f"<b>Lokasi :</b> {lokasi}", styles["Normal"]))
                        elements.append(Paragraph(f"<b>STA :</b> {sta_umum}", styles["Normal"]))
                        elements.append(Paragraph(f"<b>Surveyor :</b> {surveyor}", styles["Normal"]))
                        elements.append(Paragraph(f"<b>Tanggal :</b> {tanggal}", styles["Normal"])); elements.append(PageBreak())

                        elements.append(Paragraph("<b>1. Tabel Rekapitulasi Umum</b>", styles["Heading2"]))
                        ringkasan_table = Table([["Lokasi", lokasi], ["STA", sta_umum], ["Jumlah Segmen", str(len(seg_gdf))], ["Panjang Jalan Terukur", f"{len(seg_gdf)*interval_segmen} meter"], ["Rata-rata PCI Keseluruhan", f"{rata_pci}"], ["Kondisi Dominan", kondisi_dominan]], colWidths=[200, 400])
                        ringkasan_table.setStyle(TableStyle([('GRID',(0,0),(-1,-1),0.5,colors.grey), ('BACKGROUND',(0,0),(0,-1),colors.HexColor("#f3f4f6")), ('FONTNAME', (0,0), (0,-1), 'Helvetica-Bold'), ('PADDING', (0,0), (-1,-1), 8)]))
                        elements.append(ringkasan_table); elements.append(Spacer(1, 0.3 * inch))

                        elements.append(Paragraph("<b>2. Grafik Distribusi PCI</b>", styles["Heading2"]))
                        elements.append(Image(grafik_path, width=6.5*inch, height=3.5*inch)); elements.append(PageBreak())
                        elements.append(Paragraph("<b>3. Peta Kondisi Jalan</b>", styles["Heading2"]))
                        
                        img_reader = ImageReader(peta_path)
                        img_w, img_h = img_reader.getSize()
                        aspect = img_h / float(img_w)
                        max_w = 9.5 * inch
                        max_h = 5.5 * inch
                        
                        if max_w * aspect <= max_h:
                            final_w = max_w
                            final_h = max_w * aspect
                        else:
                            final_h = max_h
                            final_w = max_h / aspect
                            
                        elements.append(Image(peta_path, width=final_w, height=final_h))
                        elements.append(PageBreak())

                        elements.append(Paragraph("<b>LAMPIRAN: KERTAS KERJA PER SEGMEN</b>", styles["Heading1"]))
                        COLOR_HEADER_BG = colors.HexColor("#1e293b"); COLOR_HEADER_TXT = colors.white; COLOR_CELL_BG = colors.HexColor("#f8fafc")

                        for idx, seg in df_pci.sort_values('Segmen').reset_index(drop=True).iterrows():
                            if idx > 0: elements.append(PageBreak())
                            seg_id = seg["Segmen"]
                            df_seg_detail = df_detail[df_detail["Segmen"] == seg_id] if not df_detail.empty else pd.DataFrame()
                            hdv_val = df_seg_detail["DV"].max() if not df_seg_detail.empty else 0.0
                            m_val = min(1 + (9.0 / 95.0) * (100.0 - hdv_val), 10.0) if hdv_val > 0 else 0.0

                            elements.append(Paragraph(f"<b>REPORT SEGMEN : {seg_id} (STA: {seg['STA']})</b>", styles["Heading2"]))
                            elements.append(Paragraph("<b>A. Flexible Pavement Condition Data Sheet</b>", header_style))
                            tabel_a_data = [["Distress Type", "Severity", "Quantity (sq.m)", "Density (%)", "Deduct Value (DV)"]]
                            if df_seg_detail.empty: tabel_a_data.append(["Tidak ada kerusakan", "-", "0.00", "0.00", "0.00"])
                            else:
                                for _, row in df_seg_detail.iterrows():
                                    tqty = (row["Density"] / 100.0) * (interval_segmen * lebar_jalan)
                                    tabel_a_data.append([row["Distress"].replace("_", " ").title(), row["Severity"], f"{tqty:.2f}", f"{row['Density']:.2f}", f"{row['DV']:.2f}"])
                            t_a = Table(tabel_a_data, colWidths=[3*inch, 1.5*inch, 1.5*inch, 1.5*inch, 2*inch])
                            t_a.setStyle(TableStyle([('GRID', (0,0), (-1,-1), 0.5, colors.lightgrey), ('BACKGROUND', (0,0), (-1,0), COLOR_HEADER_BG), ('TEXTCOLOR', (0,0), (-1,0), COLOR_HEADER_TXT), ('ALIGN', (0,0), (-1,-1), 'CENTER'), ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold')]))
                            elements.append(t_a); elements.append(Spacer(1, 0.3*inch))

                            elements.append(Paragraph("<b>B. Maximum allowable number of distresses (m)</b>", header_style))
                            t_b = Table([["Highest Deduct Value (HDV)", "m = 1 + (9/95)*(100 - HDV) \u2264 10"], [f"{hdv_val:.2f}", f"{m_val:.2f}"]], colWidths=[4.75*inch, 4.75*inch])
                            t_b.setStyle(TableStyle([('BOX', (0,0), (-1,-1), 1, colors.lightgrey), ('INNERGRID', (0,0), (-1,-1), 0.5, colors.lightgrey), ('BACKGROUND', (0,0), (-1,0), COLOR_HEADER_BG), ('TEXTCOLOR', (0,0), (-1,0), COLOR_HEADER_TXT), ('ALIGN', (0,0), (-1,-1), 'CENTER'), ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold')]))
                            elements.append(t_b); elements.append(Spacer(1, 0.3*inch))

                            elements.append(Paragraph("<b>C. Calculate Pavement Condition Index (PCI)</b>", header_style))
                            t_c = Table([["Max CDV", "PCI = 100 - Max CDV", "Rating (ASTM)"], [f"{seg['CDV']:.2f}", f"{seg['PCI']:.2f}", seg['Rating']]], colWidths=[3.16*inch, 3.16*inch, 3.18*inch])
                            t_c.setStyle(TableStyle([('BOX', (0,0), (-1,-1), 1, colors.lightgrey), ('INNERGRID', (0,0), (-1,-1), 0.5, colors.lightgrey), ('BACKGROUND', (0,0), (-1,0), COLOR_HEADER_BG), ('TEXTCOLOR', (0,0), (-1,0), COLOR_HEADER_TXT), ('ALIGN', (0,0), (-1,-1), 'CENTER'), ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'), ('BACKGROUND', (2,1), (2,1), colors.HexColor(warna_pci.get(seg['Rating'], "#FFF")))]))
                            elements.append(t_c)
                        doc.build(elements)

                        # EXPORT GPKG & EXCEL
                        gpkg_path = os.path.join(tmpdir, "Peta_Hasil_PCI.gpkg")
                        export_gdf = seg_gdf.copy()
                        for col in export_gdf.columns:
                            if export_gdf[col].apply(lambda x: isinstance(x, (list, tuple))).any(): export_gdf[col] = export_gdf[col].astype(str)
                        export_gdf.to_file(gpkg_path, driver="GPKG")

                        excel_buffer = io.BytesIO()
                        with pd.ExcelWriter(excel_buffer, engine='xlsxwriter') as writer:
                            df_pci.to_excel(writer, sheet_name='Rekap PCI', index=False)
                            df_detail.to_excel(writer, sheet_name='Detail Kerusakan', index=False)
                        
                        # SIMPAN KE STATE
                        st.session_state.df_pci = df_pci
                        st.session_state.df_detail_pci = df_detail
                        st.session_state.seg_gdf_pci = seg_gdf           
                        st.session_state.excel_bytes_pci = excel_buffer.getvalue()   
                        with open(peta_path, "rb") as f: st.session_state.peta_bytes_pci = f.read()
                        with open(grafik_path, "rb") as f: st.session_state.grafik_bytes_pci = f.read()
                        with open(pdf_path, "rb") as f: st.session_state.pdf_bytes_pci = f.read()
                        with open(gpkg_path, "rb") as f: st.session_state.gpkg_bytes_pci = f.read()       
                        
                        st.session_state.pci_selesai = True

                    except Exception as e:
                        st.error(f"❌ Terjadi kesalahan: {e}"); st.session_state.pci_selesai = False

    # HASIL PCI TAMPILAN
    if st.session_state.pci_selesai:
        st.success("✅ Analisis PCI Berhasil!")
        col_res1, col_res2 = st.columns([2, 1])
        with col_res1:
            st.subheader("🗺️ Peta Kondisi PCI")
            
            # Siapkan base layer jalan (Aman dari NaN agar tidak crash)
            map_gdf = st.session_state.seg_gdf_pci[['geometry', 'Segmen', 'STA', 'PCI', 'Rating']].copy()
            for col in ['Segmen', 'STA', 'PCI', 'Rating']:
                map_gdf[col] = map_gdf[col].fillna("-").astype(str)
            map_gdf = map_gdf.to_crs(epsg=4326)
            
            m = folium.Map(location=[map_gdf.geometry.centroid.y.mean(), map_gdf.geometry.centroid.x.mean()], zoom_start=15)
            
            # 1. Layer Base Segmen Jalan (Kotak Presisi)
            warna_pci_dict = {"Good": "#006400", "Satisfactory": "#8FBC8F", "Fair": "#FFFF00", "Poor": "#FF6347", "Very Poor": "#FF4500", "Serious": "#8B0000", "Failed": "#A9A9A9"}
            folium.GeoJson(map_gdf, name="Segmen Jalan (PCI)",
                           style_function=lambda f: {'fillColor': warna_pci_dict.get(f['properties']['Rating'], "#000"), 'color': 'black', 'weight': 1, 'fillOpacity': 0.7},
                           tooltip=folium.features.GeoJsonTooltip(fields=['Segmen', 'STA', 'PCI', 'Rating'])).add_to(m)
            
            # 2. Layer Kerusakan Asli (Raw Polygon yang aman)
            if st.session_state.master_distress_pci is not None and not st.session_state.master_distress_pci.empty:
                # Filter kolom agar tidak crash di Folium Tooltip
                distress_map = st.session_state.master_distress_pci[['geometry', 'Distress_Type', 'Severity']].copy()
                distress_map['Distress_Type'] = distress_map['Distress_Type'].fillna("-").astype(str)
                distress_map['Severity'] = distress_map['Severity'].fillna("-").astype(str)
                distress_map = distress_map.to_crs(epsg=4326)
                
                folium.GeoJson(
                    distress_map, 
                    name="Data Kerusakan",
                    style_function=lambda x: {'color': '#e74c3c', 'fillColor': '#e74c3c', 'weight': 2, 'fillOpacity': 0.6},
                    tooltip=folium.features.GeoJsonTooltip(fields=['Distress_Type', 'Severity'])
                ).add_to(m)
                
            folium.LayerControl().add_to(m)
            st_folium(m, use_container_width=True, height=400)
        with col_res2:
            st.subheader("Distribusi")
            st.image(st.session_state.grafik_bytes_pci)
            st.metric("Rata-rata PCI", round(st.session_state.df_pci["PCI"].mean(), 2))
        
        # --- TABEL REKAP KESELURUHAN PCI ---
        st.markdown("---")
        st.subheader("📋 Tabel Rekapitulasi Keseluruhan (PCI)")
        st.caption("Klik pada header kolom (misal: 'Skor PCI') untuk mengurutkan data dari nilai tertinggi ke terendah atau sebaliknya.")
        df_pci_display = st.session_state.df_pci[['Segmen', 'STA', 'PCI', 'Rating']].copy()
        df_pci_display.columns = ['Segmen', 'STA', 'Skor PCI', 'Kelas Kerusakan']
        st.dataframe(df_pci_display, use_container_width=True, hide_index=True)
        # -----------------------------------

        # DASHBOARD PER SEGMEN (PCI)
        st.markdown("---")
        st.subheader("🔎 Dashboard Detail Per Segmen (PCI)")
        pilihan_segmen_pci = st.selectbox("Pilih Segmen:", st.session_state.df_pci["Segmen"].tolist(), key="sel_pci")

        if pilihan_segmen_pci:
            seg_data = st.session_state.df_pci[st.session_state.df_pci["Segmen"] == pilihan_segmen_pci].iloc[0]
            df_seg_detail = st.session_state.df_detail_pci[st.session_state.df_detail_pci["Segmen"] == pilihan_segmen_pci]
            hdv_val = df_seg_detail["DV"].max() if not df_seg_detail.empty else 0.0
            m_val = min(1 + (9.0 / 95.0) * (100.0 - hdv_val), 10.0) if hdv_val > 0 else 0.0

            st.markdown(f"#### REPORT SEGMEN : {pilihan_segmen_pci} (STA: {seg_data['STA']})")
            st.markdown("**A. Flexible Pavement Condition Data Sheet**")
            if df_seg_detail.empty: st.info("✅ Tidak ada kerusakan pada segmen ini.")
            else:
                display_df = df_seg_detail.copy()
                display_df["Distress Type"] = display_df["Distress"].str.replace("_", " ").str.title()
                display_df["Quantity (sq.m)"] = (display_df["Density"] / 100.0) * (interval_segmen * lebar_jalan)
                display_df = display_df[["Distress Type", "Severity", "Quantity (sq.m)", "Density", "DV"]].rename(columns={"Density": "Density (%)", "DV": "Deduct Value (DV)"})
                st.dataframe(display_df.style.format({"Quantity (sq.m)": "{:.2f}", "Density (%)": "{:.2f}", "Deduct Value (DV)": "{:.2f}"}), use_container_width=True, hide_index=True)

            st.markdown("**B. Maximum allowable number of distresses (m)**")
            col_b1, col_b2 = st.columns(2)
            with col_b1: st.markdown(metric_card("Highest Deduct Value (HDV)", f"{hdv_val:.2f}"), unsafe_allow_html=True)
            with col_b2: st.markdown(metric_card("m = 1 + (9/95)*(100 - HDV) ≤ 10", f"{m_val:.2f}"), unsafe_allow_html=True)

            st.markdown("<br>**C. Calculate Pavement Condition Index (PCI)**", unsafe_allow_html=True)
            col_c1, col_c2, col_c3 = st.columns(3)
            with col_c1: st.markdown(metric_card("Max CDV", f"{seg_data['CDV']:.2f}"), unsafe_allow_html=True)
            with col_c2: st.markdown(metric_card("PCI = 100 - Max CDV", f"{seg_data['PCI']:.2f}"), unsafe_allow_html=True)
            bg_col = warna_pci_dict.get(seg_data['Rating'], "#FFFFFF")
            txt_col = "#000000" if seg_data['Rating'] in ["Satisfactory", "Fair", "Good"] else "#ffffff"
            with col_c3: st.markdown(metric_card("Rating (ASTM)", seg_data['Rating'], value_color=txt_col, bg_color=bg_col, text_color=txt_col), unsafe_allow_html=True)

        st.markdown("---")
        st.subheader("💾 Download Hasil Analisis PCI")
        col_dl1, col_dl2, col_dl3 = st.columns(3)
        with col_dl1: st.download_button("📄 Laporan Full PDF", data=st.session_state.pdf_bytes_pci, file_name=f"Laporan_PCI_{lokasi.replace(' ', '_')}.pdf", mime="application/pdf", type="primary", use_container_width=True)
        with col_dl2: st.download_button("🗺️ Peta Spasial (.gpkg)", data=st.session_state.gpkg_bytes_pci, file_name=f"Peta_PCI_{lokasi.replace(' ', '_')}.gpkg", mime="application/geopackage+sqlite3", type="secondary", use_container_width=True)
        with col_dl3: st.download_button("📊 Data Mentah (.xlsx)", data=st.session_state.excel_bytes_pci, file_name=f"Data_PCI_{lokasi.replace(' ', '_')}.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", type="secondary", use_container_width=True)

# =========================================
# MODUL SDI
# =========================================
elif menu == "📉 Modul SDI (Surface Distress Index)":
    st.title("🛣️ Modul Evaluasi SDI")
    st.markdown("Otomatisasi perhitungan Surface Distress Index (SDI) berbasis GIS.")
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.subheader("📁 1. Data Jalan")
        jalan_file_sdi = st.file_uploader("Upload Shapefile Jalan (.zip)", type="zip", key="sdi_jalan")
    with col2:
        st.subheader("⚠️ 2. Data Kerusakan")
        retak_file = st.file_uploader("SHP Retak (.zip)", type="zip", key="sdi_retak")
        pothole_file = st.file_uploader("SHP Lubang (.zip)", type="zip", key="sdi_pothole")
        rutting_file = st.file_uploader("SHP Rutting (.zip)", type="zip", key="sdi_rutting")
    with col3:
        st.subheader("🗺️ 3. DSM")
        dsm_mode_sdi = st.radio("Cara Input Data DSM:", ["Upload File .tif", "Paste Link Google Drive"], key="sdi_dsm_mode")
        dsm_file_sdi = None; dsm_link_sdi = ""
        if dsm_mode_sdi == "Upload File .tif": dsm_file_sdi = st.file_uploader("Upload Data DSM (.tif)", type="tif", key="sdi_dsm_file")
        else: dsm_link_sdi = st.text_input("Paste Link Drive (.tif)", key="sdi_dsm_link")

    if st.button("🚀 Proses & Hitung SDI", type="primary", use_container_width=True):
        is_dsm_valid = (dsm_mode_sdi == "Upload File .tif" and dsm_file_sdi is not None) or (dsm_mode_sdi == "Paste Link Google Drive" and dsm_link_sdi != "")
        if not jalan_file_sdi or not is_dsm_valid: st.error("⚠️ Mohon lengkapi Shapefile Jalan dan Data DSM.")
        else:
            with st.spinner("Memproses Analisis SDI..."):
                with tempfile.TemporaryDirectory() as tmpdir:
                    try:
                        dsm_path = os.path.join(tmpdir, "dsm.tif")
                        if dsm_mode_sdi == "Upload File .tif":
                            with open(dsm_path, "wb") as f: f.write(dsm_file_sdi.getbuffer())
                        else:
                            import gdown, re
                            match = re.search(r"/d/([a-zA-Z0-9_-]+)", dsm_link_sdi)
                            if match: gdown.download(id=match.group(1), output=dsm_path, quiet=False)

                        jalan = read_zip_shapefile(jalan_file_sdi, tmpdir)
                        if jalan.crs is None: jalan.set_crs(epsg=4326, inplace=True)
                        if jalan.crs.to_epsg() != epsg_code: jalan = jalan.to_crs(epsg=epsg_code)
                        union_geom = jalan.geometry.union_all()
                        merged_line = linemerge(union_geom) if union_geom.geom_type == "MultiLineString" else union_geom
                        panjang_total = merged_line.length
                        segments = [substring(merged_line, start, min(start + interval_segmen, panjang_total)) 
                                    for start in np.arange(0, panjang_total, interval_segmen)]
                        seg_gdf = gpd.GeoDataFrame(geometry=segments, crs=jalan.crs)
                        seg_gdf["Segmen"] = range(1, len(seg_gdf)+1)
                        seg_gdf["STA"] = seg_gdf["Segmen"].apply(lambda x: f"{(x-1)*interval_segmen:03.0f}+000 - {min(x*interval_segmen, int(panjang_total)):03.0f}+000")
                        seg_gdf["geometry"] = seg_gdf.buffer(lebar_jalan / 2, cap_style=2)
                        seg_gdf["Luas_Segmen"] = seg_gdf.geometry.area
                        
                        gdf_retak = read_zip_shapefile(retak_file, tmpdir) if retak_file else gpd.GeoDataFrame(columns=['geometry'], crs=seg_gdf.crs)
                        gdf_pothole = read_zip_shapefile(pothole_file, tmpdir) if pothole_file else gpd.GeoDataFrame(columns=['geometry'], crs=seg_gdf.crs)
                        gdf_rutting = read_zip_shapefile(rutting_file, tmpdir) if rutting_file else gpd.GeoDataFrame(columns=['geometry'], crs=seg_gdf.crs)
                        
                        for gdf in [gdf_retak, gdf_pothole, gdf_rutting]:
                            if not gdf.empty:
                                if gdf.crs is None: gdf.set_crs(seg_gdf.crs, inplace=True)
                                elif gdf.crs != seg_gdf.crs: gdf.to_crs(seg_gdf.crs, inplace=True)

                        if not gdf_rutting.empty: gdf_rutting = hitung_depth_cm(gdf_rutting, dsm_path)

                        hasil_sdi = []
                        for idx, seg in seg_gdf.iterrows():
                            seg_poly = gpd.GeoDataFrame(geometry=[seg.geometry], crs=seg_gdf.crs)
                            luas_seg = seg["Luas_Segmen"]
                            
                            persen_retak = 0.0; lebar_retak = 0.0
                            if not gdf_retak.empty:
                                retak_seg = gpd.overlay(gdf_retak, seg_poly, how="intersection")
                                if not retak_seg.empty:
                                    luas_retak = retak_seg.geometry.area.sum()
                                    persen_retak = (luas_retak / luas_seg) * 100 if luas_seg > 0 else 0
                                    lengths = retak_seg.geometry.length
                                    valid_lengths = lengths[lengths > 0]
                                    if len(valid_lengths) > 0:
                                        retak_seg.loc[lengths > 0, "lebar_calc"] = retak_seg.geometry.area / valid_lengths
                                        lebar_retak = retak_seg["lebar_calc"].mean() * 1000 

                            jumlah_lubang = 0
                            if not gdf_pothole.empty:
                                pothole_seg = gpd.sjoin(gdf_pothole, seg_poly, predicate="within")
                                jumlah_lubang = len(pothole_seg)

                            kedalaman_rutting = 0.0
                            if not gdf_rutting.empty:
                                rutting_seg = gpd.overlay(gdf_rutting, seg_poly, how="intersection")
                                if not rutting_seg.empty: kedalaman_rutting = rutting_seg["kedalaman_calc"].mean()

                            kedalaman_rutting = 0 if pd.isna(kedalaman_rutting) else kedalaman_rutting
                            sdi1, sdi2, sdi3, sdi4, kondisi = hitung_sdi(persen_retak, lebar_retak, jumlah_lubang, kedalaman_rutting)
                            
                            hasil_sdi.append({"Segmen": seg["Segmen"], "STA": seg["STA"], "%Retak": round(persen_retak, 2), "Lebar Retak (mm)": round(lebar_retak, 2), "Jumlah Lubang": jumlah_lubang, "Rutting (cm)": round(kedalaman_rutting, 2), "SDI1": sdi1, "SDI2": sdi2, "SDI3": sdi3, "SDI": round(sdi4, 2), "Kondisi": kondisi})

                        df_sdi = pd.DataFrame(hasil_sdi)
                        seg_gdf = seg_gdf.merge(df_sdi.drop(columns=["STA"]), on="Segmen", how="left")
                        
                        warna_kondisi = {"Baik": "#2ecc71", "Sedang": "#f1c40f", "Rusak Ringan": "#e67e22", "Rusak Berat": "#e74c3c"}
                        fig_map, ax_map = plt.subplots(figsize=(10,6))
                        seg_gdf.boundary.plot(ax=ax_map, linewidth=0.5, color="black")
                        legend_handles = []
                        for kondisi, warna in warna_kondisi.items():
                            subset = seg_gdf[seg_gdf["Kondisi"] == kondisi]
                            if not subset.empty: 
                                subset.plot(ax=ax_map, color=warna, edgecolor="black", linewidth=1)
                                legend_handles.append(mpatches.Patch(color=warna, label=f"{kondisi} ({len(subset)})"))
                        for idx, row in seg_gdf.iterrows():
                            centroid = row.geometry.centroid
                            ax_map.text(centroid.x, centroid.y, f"S{row['Segmen']}\n{row['SDI']:.0f}", fontsize=7, weight="bold", ha="center", va="center", bbox=dict(facecolor="white", alpha=0.8, boxstyle="round,pad=0.2", edgecolor="gray", lw=0.5))
                        if legend_handles: ax_map.legend(handles=legend_handles, loc="best", title="Kategori Kondisi", fontsize=8)
                        ax_map.axis("off")
                        peta_path = os.path.join(tmpdir, "peta_sdi.png")
                        plt.savefig(peta_path, dpi=300, bbox_inches='tight')
                        plt.close(fig_map)
                        
                        fig_bar, ax_bar = plt.subplots(figsize=(6,4))
                        rekap = seg_gdf["Kondisi"].value_counts()
                        rekap.plot(kind="bar", color=[warna_kondisi.get(x, "grey") for x in rekap.index], edgecolor="black", ax=ax_bar)
                        plt.title("Distribusi Kondisi Jalan (Segmen)")
                        plt.xticks(rotation=0)
                        plt.tight_layout()
                        grafik_path = os.path.join(tmpdir, "grafik_sdi.png")
                        plt.savefig(grafik_path, dpi=300)
                        plt.close(fig_bar)
                        
                        # Generate Full PDF Report (SDI)
                        pdf_path = os.path.join(tmpdir, "Laporan_SDI.pdf")
                        doc = SimpleDocTemplate(pdf_path, pagesize=pagesizes.A4, rightMargin=30, leftMargin=30, topMargin=30, bottomMargin=30)
                        elements = []
                        styles = getSampleStyleSheet()
                        cover_style = ParagraphStyle('cover', parent=styles['Title'], alignment=TA_CENTER)
                        rata_sdi = round(df_sdi["SDI"].mean(), 2)
                        kondisi_dominan = df_sdi["Kondisi"].value_counts().idxmax() if not df_sdi.empty else "-"

                        elements.append(Paragraph(instansi, cover_style)); elements.append(Spacer(1, 0.3*inch))
                        elements.append(Paragraph("LAPORAN SURVEY", cover_style)); elements.append(Spacer(1, 0.3*inch))
                        elements.append(Paragraph("SURFACE DISTRESS INDEX (SDI)", cover_style)); elements.append(Spacer(1, 1*inch))
                        elements.append(Paragraph(f"<b>Lokasi :</b> {lokasi}", styles["Normal"]))
                        elements.append(Paragraph(f"<b>STA :</b> {sta_umum}", styles["Normal"]))
                        elements.append(Paragraph(f"<b>Surveyor :</b> {surveyor}", styles["Normal"]))
                        elements.append(Paragraph(f"<b>Tanggal :</b> {tanggal}", styles["Normal"])); elements.append(PageBreak())

                        elements.append(Paragraph("<b>1. Ringkasan Rekapitulasi Umum</b>", styles["Heading2"]))
                        ringkasan_table = Table([["Lokasi", lokasi], ["STA", sta_umum], ["Jumlah Segmen", str(len(seg_gdf))], ["Panjang Jalan Terukur", f"{len(seg_gdf)*interval_segmen} meter"], ["Rata-rata SDI Keseluruhan", f"{rata_sdi}"], ["Kondisi Dominan", kondisi_dominan]], colWidths=[200, 300])
                        ringkasan_table.setStyle(TableStyle([('GRID',(0,0),(-1,-1),0.5,colors.grey), ('BACKGROUND',(0,0),(0,-1),colors.HexColor("#f3f4f6")), ('FONTNAME', (0,0), (0,-1), 'Helvetica-Bold')]))
                        elements.append(ringkasan_table); elements.append(Spacer(1, 0.3 * inch))

                        elements.append(Paragraph("<b>2. Visualisasi Kondisi Jalan</b>", styles["Heading2"]))
                        
                        img_reader_sdi = ImageReader(peta_path)
                        img_w_sdi, img_h_sdi = img_reader_sdi.getSize()
                        aspect_sdi = img_h_sdi / float(img_w_sdi)
                        max_w_sdi = 7.5 * inch
                        max_h_sdi = 4.5 * inch
                        
                        if max_w_sdi * aspect_sdi <= max_h_sdi:
                            final_w_sdi = max_w_sdi
                            final_h_sdi = max_w_sdi * aspect_sdi
                        else:
                            final_h_sdi = max_h_sdi
                            final_w_sdi = max_h_sdi / aspect_sdi
                            
                        elements.append(Image(peta_path, width=final_w_sdi, height=final_h_sdi))
                        elements.append(Spacer(1, 0.2 * inch))
                        elements.append(Image(grafik_path, width=4.5*inch, height=3*inch))
                        elements.append(PageBreak())

                        elements.append(Paragraph("<b>3. Data Kerusakan Terukur Per Segmen</b>", styles["Heading2"]))
                        tabel1_data = [["Segmen", "STA", "% Retak", "Lebar Retak\n(mm)", "Jumlah\nLubang", "Rutting\n(cm)"]]
                        for _, row in df_sdi.iterrows():
                            tabel1_data.append([str(row["Segmen"]), row["STA"], str(row["%Retak"]), str(row["Lebar Retak (mm)"]), str(row["Jumlah Lubang"]), str(row["Rutting (cm)"])])
                        t1_detail = Table(tabel1_data, repeatRows=1, colWidths=[0.8*inch, 2.0*inch, 1.0*inch, 1.2*inch, 1.0*inch, 1.0*inch])
                        t1_detail.setStyle(TableStyle([('GRID', (0,0), (-1,-1), 0.5, colors.grey), ('BACKGROUND', (0,0), (-1,0), colors.HexColor("#1e293b")), ('TEXTCOLOR', (0,0), (-1,0), colors.white), ('ALIGN', (0,0), (-1,-1), 'CENTER'), ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'), ('FONTSIZE', (0,0), (-1,-1), 9)]))
                        elements.append(t1_detail); elements.append(PageBreak())

                        elements.append(Paragraph("<b>4. Perhitungan Berjenjang SDI Per Segmen</b>", styles["Heading2"]))
                        tabel2_data = [["Segmen", "STA", "SDI 1\n(Retak)", "SDI 2\n(+L. Retak)", "SDI 3\n(+Lubang)", "SDI\n(+Rutting)", "Kondisi Akhir"]]
                        for _, row in df_sdi.iterrows():
                            tabel2_data.append([str(row["Segmen"]), row["STA"], str(row["SDI1"]), str(row["SDI2"]), str(row["SDI3"]), str(row["SDI"]), row["Kondisi"]])
                        t2_detail = Table(tabel2_data, repeatRows=1, colWidths=[0.8*inch, 1.8*inch, 0.8*inch, 0.9*inch, 0.8*inch, 0.8*inch, 1.1*inch])
                        t2_detail.setStyle(TableStyle([('GRID', (0,0), (-1,-1), 0.5, colors.grey), ('BACKGROUND', (0,0), (-1,0), colors.HexColor("#1e293b")), ('TEXTCOLOR', (0,0), (-1,0), colors.white), ('ALIGN', (0,0), (-1,-1), 'CENTER'), ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'), ('FONTSIZE', (0,0), (-1,-1), 9)]))
                        elements.append(t2_detail)
                        doc.build(elements)

                        gpkg_path = os.path.join(tmpdir, "Peta_Hasil_SDI.gpkg")
                        export_gdf = seg_gdf.copy()
                        for col in export_gdf.columns:
                            if export_gdf[col].apply(lambda x: isinstance(x, (list, tuple))).any(): export_gdf[col] = export_gdf[col].astype(str)
                        export_gdf.to_file(gpkg_path, driver="GPKG")

                        excel_buffer = io.BytesIO()
                        with pd.ExcelWriter(excel_buffer, engine='xlsxwriter') as writer: df_sdi.to_excel(writer, sheet_name='Rekap SDI', index=False)
                        
                        st.session_state.df_sdi = df_sdi
                        st.session_state.seg_gdf_sdi = seg_gdf           
                        st.session_state.excel_bytes_sdi = excel_buffer.getvalue()   
                        with open(peta_path, "rb") as f: st.session_state.peta_bytes_sdi = f.read()
                        with open(grafik_path, "rb") as f: st.session_state.grafik_bytes_sdi = f.read()
                        with open(pdf_path, "rb") as f: st.session_state.pdf_bytes_sdi = f.read()
                        with open(gpkg_path, "rb") as f: st.session_state.gpkg_bytes_sdi = f.read()       
                            
                        st.session_state.gdf_retak_sdi = gdf_retak.copy()
                        st.session_state.gdf_pothole_sdi = gdf_pothole.copy()
                        st.session_state.gdf_rutting_sdi = gdf_rutting.copy()
                        
                        st.session_state.sdi_selesai = True

                    except Exception as e: st.error(f"❌ Terjadi kesalahan: {e}"); st.session_state.sdi_selesai = False

    if st.session_state.sdi_selesai:
        st.success("✅ Analisis SDI Berhasil!")
        col_res1, col_res2 = st.columns([2, 1])
        with col_res1:
            st.subheader("🗺️ Peta Kondisi SDI")
            
            # Siapkan base layer jalan (Aman dari NaN)
            map_gdf = st.session_state.seg_gdf_sdi[['geometry', 'Segmen', 'STA', 'SDI', 'Kondisi']].copy()
            for col in ['Segmen', 'STA', 'SDI', 'Kondisi']:
                map_gdf[col] = map_gdf[col].fillna("-").astype(str)
            map_gdf = map_gdf.to_crs(epsg=4326)
            
            m = folium.Map(location=[map_gdf.geometry.centroid.y.mean(), map_gdf.geometry.centroid.x.mean()], zoom_start=15)
            
            # 1. Layer Base Segmen Jalan SDI (Kotak Presisi)
            warna_kondisi_dict = {"Baik": "#2ecc71", "Sedang": "#f1c40f", "Rusak Ringan": "#e67e22", "Rusak Berat": "#e74c3c"}
            folium.GeoJson(map_gdf, name="Segmen Jalan (SDI)",
                           style_function=lambda f: {'fillColor': warna_kondisi_dict.get(f['properties']['Kondisi'], "#000"), 'color': 'black', 'weight': 1, 'fillOpacity': 0.7},
                           tooltip=folium.features.GeoJsonTooltip(fields=['Segmen', 'STA', 'SDI', 'Kondisi'])).add_to(m)
            
            # 2. Layer Overlay Retak (Hanya ambil geometry agar tidak memicu error)
            if st.session_state.gdf_retak_sdi is not None and not st.session_state.gdf_retak_sdi.empty:
                clean_retak = st.session_state.gdf_retak_sdi[['geometry']].to_crs(epsg=4326)
                folium.GeoJson(clean_retak, name="Retak", style_function=lambda x: {'color': '#e74c3c', 'weight': 2}).add_to(m)
            
            # 3. Layer Overlay Lubang
            if st.session_state.gdf_pothole_sdi is not None and not st.session_state.gdf_pothole_sdi.empty:
                clean_pothole = st.session_state.gdf_pothole_sdi[['geometry']].to_crs(epsg=4326)
                folium.GeoJson(clean_pothole, name="Lubang", style_function=lambda x: {'color': '#3498db', 'fillColor': '#3498db', 'weight': 2, 'fillOpacity': 0.6}).add_to(m)
                               
            # 4. Layer Overlay Rutting
            if st.session_state.gdf_rutting_sdi is not None and not st.session_state.gdf_rutting_sdi.empty:
                clean_rutting = st.session_state.gdf_rutting_sdi[['geometry']].to_crs(epsg=4326)
                folium.GeoJson(clean_rutting, name="Rutting", style_function=lambda x: {'color': '#9b59b6', 'fillColor': '#9b59b6', 'weight': 2, 'fillOpacity': 0.6}).add_to(m)

            folium.LayerControl().add_to(m)
            st_folium(m, use_container_width=True, height=400)
        with col_res2:
            st.subheader("Distribusi")
            st.image(st.session_state.grafik_bytes_sdi)
            st.metric("Rata-rata Nilai SDI", round(st.session_state.df_sdi["SDI"].mean(), 2))
        
        # --- TABEL REKAP KESELURUHAN SDI ---
        st.markdown("---")
        st.subheader("📋 Tabel Rekapitulasi Keseluruhan (SDI)")
        st.caption("Klik pada header kolom (misal: 'Skor SDI') untuk mengurutkan data dari nilai tertinggi ke terendah atau sebaliknya.")
        df_sdi_display = st.session_state.df_sdi[['Segmen', 'STA', 'SDI', 'Kondisi']].copy()
        df_sdi_display.columns = ['Segmen', 'STA', 'Skor SDI', 'Kelas Kerusakan']
        st.dataframe(df_sdi_display, use_container_width=True, hide_index=True)
        # -----------------------------------

        # DASHBOARD PER SEGMEN (SDI)
        st.markdown("---")
        st.subheader("🔎 Dashboard Detail Perhitungan Segmen (SDI)")
        pilihan_segmen_sdi = st.selectbox("Pilih Segmen:", st.session_state.df_sdi["Segmen"].tolist(), key="sel_sdi")

        if pilihan_segmen_sdi:
            seg_data = st.session_state.df_sdi[st.session_state.df_sdi["Segmen"] == pilihan_segmen_sdi].iloc[0]

            st.markdown(f"#### REPORT SEGMEN : {pilihan_segmen_sdi} (STA: {seg_data['STA']})")
            st.markdown("**A. Data Kerusakan Terukur**")
            col_m1, col_m2, col_m3, col_m4 = st.columns(4)
            with col_m1: st.markdown(metric_card("Luas Retak (%)", f"{seg_data['%Retak']:.2f}%"), unsafe_allow_html=True)
            with col_m2: st.markdown(metric_card("Lebar Retak (mm)", f"{seg_data['Lebar Retak (mm)']:.2f}"), unsafe_allow_html=True)
            with col_m3: st.markdown(metric_card("Jumlah Lubang (Ttk)", f"{seg_data['Jumlah Lubang']}"), unsafe_allow_html=True)
            with col_m4: st.markdown(metric_card("Rutting/Alur (cm)", f"{seg_data['Rutting (cm)']:.2f}"), unsafe_allow_html=True)

            st.markdown("<br>**B. Perhitungan Berjenjang SDI**", unsafe_allow_html=True)
            col_s1, col_s2, col_s3, col_s4, col_s5 = st.columns(5)
            with col_s1: st.markdown(metric_card("SDI 1<br>(Retak)", f"{seg_data['SDI1']}"), unsafe_allow_html=True)
            with col_s2: st.markdown(metric_card("SDI 2<br>(+Lebar Retak)", f"{seg_data['SDI2']}"), unsafe_allow_html=True)
            with col_s3: st.markdown(metric_card("SDI 3<br>(+Lubang)", f"{seg_data['SDI3']}"), unsafe_allow_html=True)
            with col_s4: st.markdown(metric_card("SDI<br>(Nilai Akhir)", f"{seg_data['SDI']:.2f}", value_color="#ffcc00"), unsafe_allow_html=True)
            
            bg_col = warna_kondisi_dict.get(seg_data['Kondisi'], "#FFFFFF")
            txt_col = "#000000" if seg_data['Kondisi'] in ["Sedang", "Baik"] else "#ffffff"
            with col_s5: st.markdown(metric_card("Kondisi<br>Akhir", seg_data['Kondisi'], value_color=txt_col, bg_color=bg_col, text_color=txt_col), unsafe_allow_html=True)

        st.markdown("---")
        st.subheader("💾 Download Hasil Analisis SDI")
        col_dl1, col_dl2, col_dl3 = st.columns(3)
        with col_dl1: st.download_button("📄 Laporan Full PDF", data=st.session_state.pdf_bytes_sdi, file_name=f"Laporan_SDI_{lokasi.replace(' ', '_')}.pdf", mime="application/pdf", type="primary", use_container_width=True)
        with col_dl2: st.download_button("🗺️ Peta Spasial (.gpkg)", data=st.session_state.gpkg_bytes_sdi, file_name=f"Peta_SDI_{lokasi.replace(' ', '_')}.gpkg", mime="application/geopackage+sqlite3", type="secondary", use_container_width=True)
        with col_dl3: st.download_button("📊 Data Mentah (.xlsx)", data=st.session_state.excel_bytes_sdi, file_name=f"Data_SDI_{lokasi.replace(' ', '_')}.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", type="secondary", use_container_width=True)

# =========================================
# MODUL KOMPARASI
# =========================================
elif menu == "📊 Komparasi (PCI vs SDI)":
    st.title("📊 Dashboard Komparasi (PCI vs SDI)")
    
    if st.session_state.pci_selesai and st.session_state.sdi_selesai:
        st.success("✅ Data PCI dan SDI berhasil disinkronisasi.")
        
        # 1. PERSIAPAN DATA MERGE
        df_komparasi = pd.merge(
            st.session_state.df_pci[['Segmen', 'STA', 'PCI', 'Rating']],
            st.session_state.df_sdi[['Segmen', 'SDI', 'Kondisi']],
            on='Segmen',
            how='inner'
        )
        
        # Gabungkan dengan geometri untuk keperluan Peta
        gdf_komparasi = st.session_state.seg_gdf_pci[['Segmen', 'geometry']].merge(df_komparasi, on='Segmen')

        # ==========================================
        # 2. KARTU METRIK RINGKASAN
        # ==========================================
        st.markdown("### 📌 Ringkasan Eksekutif")
        col_m1, col_m2, col_m3 = st.columns(3)
        
        with col_m1:
            st.metric(label="Rata-rata Nilai PCI", value=round(df_komparasi["PCI"].mean(), 2))
        with col_m2:
            st.metric(label="Rata-rata Nilai SDI", value=round(df_komparasi["SDI"].mean(), 2))
        with col_m3:
            # Definisi Kritis: PCI di bawah 40 (Poor ke bawah) ATAU SDI di atas 100 (Rusak Ringan/Berat)
            segmen_kritis = len(df_komparasi[(df_komparasi["PCI"] <= 40) | (df_komparasi["SDI"] > 100)])
            st.metric(label="🚨 Segmen Kritis (Prioritas)", value=f"{segmen_kritis} Segmen")
            
        st.caption("**Kriteria Segmen Kritis:** Ditentukan apabila suatu segmen memiliki nilai **PCI ≤ 40** (Kategori *Poor, Very Poor, Serious, Failed*) **ATAU** nilai **SDI > 100** (Kategori *Rusak Ringan, Rusak Berat*).")

        st.markdown("---")

        # ==========================================
        # 3. LEGENDA RATING (URUTAN KLASIFIKASI)
        # ==========================================
        st.markdown("### 📖 Legenda Klasifikasi Kondisi Jalan")
        col_leg1, col_leg2 = st.columns(2)
        
        with col_leg1:
            st.markdown("**Urutan Pavement Condition Index (PCI)**")
            skala_pci = [
                ("Good", "#006400", "white", "86 - 100"),
                ("Satisfactory", "#8FBC8F", "black", "71 - 85"),
                ("Fair", "#FFFF00", "black", "56 - 70"),
                ("Poor", "#FF6347", "white", "41 - 55"),
                ("Very Poor", "#FF4500", "white", "26 - 40"),
                ("Serious", "#8B0000", "white", "11 - 25"),
                ("Failed", "#A9A9A9", "black", "0 - 10")
            ]
            for nama, bg, txt, rentang in skala_pci:
                st.markdown(f"<div style='background-color: {bg}; color: {txt}; padding: 5px 10px; margin-bottom: 2px; border-radius: 3px; display: flex; justify-content: space-between; font-size: 14px; font-weight: bold;'><span>{nama}</span><span>{rentang}</span></div>", unsafe_allow_html=True)
        
        with col_leg2:
            st.markdown("**Urutan Surface Distress Index (SDI)**")
            skala_sdi = [
                ("Baik", "#2ecc71", "white", "< 50"),
                ("Sedang", "#f1c40f", "black", "50 - 100"),
                ("Rusak Ringan", "#e67e22", "white", "101 - 150"),
                ("Rusak Berat", "#e74c3c", "white", "> 150")
            ]
            for nama, bg, txt, rentang in skala_sdi:
                st.markdown(f"<div style='background-color: {bg}; color: {txt}; padding: 5px 10px; margin-bottom: 2px; border-radius: 3px; display: flex; justify-content: space-between; font-size: 14px; font-weight: bold;'><span>{nama}</span><span>{rentang}</span></div>", unsafe_allow_html=True)

        st.markdown("---")

        # ==========================================
        # 4. PETA KOMPARASI (SIDE-BY-SIDE)
        # ==========================================
        st.markdown("### 🗺️ Komparasi Spasial")
        st.caption("Geser dan zoom peta di bawah ini untuk membandingkan sebaran kerusakan secara visual.")
        
        col_map1, col_map2 = st.columns(2)

        map_gdf = gdf_komparasi.to_crs(epsg=4326)
        center_y = map_gdf.geometry.centroid.y.mean()
        center_x = map_gdf.geometry.centroid.x.mean()

        # Peta Kiri: PCI
        with col_map1:
            st.markdown("**1. Peta Pavement Condition Index (PCI)**")
            m_pci = folium.Map(location=[center_y, center_x], zoom_start=15, tiles="CartoDB positron")
            warna_pci_dict = {"Good": "#006400", "Satisfactory": "#8FBC8F", "Fair": "#FFFF00", "Poor": "#FF6347", "Very Poor": "#FF4500", "Serious": "#8B0000", "Failed": "#A9A9A9"}
            folium.GeoJson(
                map_gdf, 
                style_function=lambda f: {'fillColor': warna_pci_dict.get(f['properties']['Rating'], "#000"), 'color': 'black', 'weight': 1, 'fillOpacity': 0.8},
                tooltip=folium.features.GeoJsonTooltip(fields=['Segmen', 'PCI', 'Rating'])
            ).add_to(m_pci)
            st_folium(m_pci, use_container_width=True, height=400, key="map_comp_pci")

        # Peta Kanan: SDI
        with col_map2:
            st.markdown("**2. Peta Surface Distress Index (SDI)**")
            m_sdi = folium.Map(location=[center_y, center_x], zoom_start=15, tiles="CartoDB positron")
            warna_kondisi_dict = {"Baik": "#2ecc71", "Sedang": "#f1c40f", "Rusak Ringan": "#e67e22", "Rusak Berat": "#e74c3c"}
            folium.GeoJson(
                map_gdf, 
                style_function=lambda f: {'fillColor': warna_kondisi_dict.get(f['properties']['Kondisi'], "#000"), 'color': 'black', 'weight': 1, 'fillOpacity': 0.8},
                tooltip=folium.features.GeoJsonTooltip(fields=['Segmen', 'SDI', 'Kondisi'])
            ).add_to(m_sdi)
            st_folium(m_sdi, use_container_width=True, height=400, key="map_comp_sdi")

        st.markdown("---")

        # ==========================================
        # 5. GRAFIK KORELASI & ANALISIS KUADRAN
        # ==========================================
        st.markdown("### 📈 Analisis Lanjutan & Distribusi Kuadran")
        st.caption("Menganalisis tren hubungan antara nilai struktural (PCI) dan kerusakan permukaan (SDI).")
        
        # Fungsi Penentuan Kuadran Otomatis
        def quadrant_label(pci, sdi):
            if pci < 55 and sdi > 100: return 'Kritis (Prioritas Utama)'
            if pci >= 55 and sdi > 100: return 'Permukaan Buruk'
            if pci < 55 and sdi <= 100: return 'Struktur Buruk'
            return 'Kondisi Baik'
            
        # Tambahkan kolom kuadran ke dataframe
        df_komparasi['Kuadran'] = df_komparasi.apply(lambda r: quadrant_label(r.PCI, r.SDI), axis=1)
        
        col_chart1, col_chart2 = st.columns(2)
        
        # --- BAGIAN KIRI: SCATTER PLOT & REGRESI ---
        with col_chart1:
            st.markdown("**Scatter Plot (PCI vs SDI)**")
            fig, ax = plt.subplots(figsize=(6, 5))
            fig.patch.set_facecolor('none')
            ax.set_facecolor('none')
            
            # Scatter titik data
            ax.scatter(df_komparasi["PCI"], df_komparasi["SDI"], color='#00a4d6', alpha=0.9, edgecolors='white', s=80, zorder=3)
            
            # Hitung Regresi & Pearson Correlation (Jika data lebih dari 1)
            if len(df_komparasi) > 1:
                from scipy.stats import pearsonr
                slope, intercept = np.polyfit(df_komparasi['PCI'], df_komparasi['SDI'], 1)
                r_val, p_val = pearsonr(df_komparasi['PCI'], df_komparasi['SDI'])
                r_squared = r_val**2
                
                # Plot Garis Tren
                x_vals = np.array(ax.get_xlim())
                y_vals = intercept + slope * x_vals
                ax.plot(x_vals, y_vals, '--', color='#2ecc71', alpha=0.8, label=f'Trend (R²={r_squared:.2f})')
                
                # Tampilkan KPI Statistik di Streamlit
                st.info(f"**Statistik:** Pearson r = `{r_val:.3f}` | p-value = `{p_val:.3e}` | R² = `{r_squared:.3f}`")
            
            # Garis batas Kritis
            ax.axvline(x=55, color='#e74c3c', linestyle='-', alpha=0.5, label='Batas Kritis PCI (55)')
            ax.axhline(y=100, color='#f39c12', linestyle='-', alpha=0.5, label='Batas Kritis SDI (100)')
            
            # Format Tema Gelap
            ax.set_xlabel("Nilai PCI (↑ Semakin Baik)", color='white')
            ax.set_ylabel("Nilai SDI (↓ Semakin Baik)", color='white')
            ax.tick_params(colors='white')
            for spine in ax.spines.values(): spine.set_edgecolor('lightgray')
            ax.grid(True, linestyle=':', alpha=0.3, color='lightgray', zorder=0)
            
            legend = ax.legend(loc="upper right", fontsize=8, facecolor='#1e293b', edgecolor='gray')
            for text in legend.get_texts(): text.set_color("white")
            
            st.pyplot(fig)

        # --- BAGIAN KANAN: DONUT CHART KUADRAN ---
        with col_chart2:
            st.markdown("**Distribusi Kuadran Kerusakan**")
            fig_pie, ax_pie = plt.subplots(figsize=(6, 5))
            fig_pie.patch.set_facecolor('none')
            
            kuadran_counts = df_komparasi['Kuadran'].value_counts()
            
            # Mapping warna sesuai kuadran
            warna_kuadran = {
                'Kritis (Prioritas Utama)': '#e74c3c',  # Merah
                'Permukaan Buruk': '#f39c12',           # Oranye
                'Struktur Buruk': '#9b59b6',            # Ungu
                'Kondisi Baik': '#2ecc71'               # Hijau
            }
            colors_pie = [warna_kuadran.get(x, '#gray') for x in kuadran_counts.index]
            
            # Buat parameter Donut Chart
            wedges, texts, autotexts = ax_pie.pie(
                kuadran_counts, labels=kuadran_counts.index, autopct='%1.1f%%',
                colors=colors_pie, startangle=90, textprops={'color': "white"},
                wedgeprops=dict(width=0.4, edgecolor='white') # Mengatur lebar cincin donut
            )
            
            plt.setp(autotexts, size=10, weight="bold")
            ax_pie.axis('equal') # Pastikan bentuknya bulat sempurna
            st.pyplot(fig_pie)

        st.markdown("---")
        st.markdown("### 📋 Tabel Detail Komparasi")
        st.caption("Data mentah hasil penggabungan analisis PCI dan SDI.")
        # Menghapus sementara kolom Kuadran dari tampilan tabel agar tidak terlalu panjang
        st.dataframe(df_komparasi.drop(columns=['Kuadran']), use_container_width=True, hide_index=True)

    else:
        st.warning("⚠️ Data belum lengkap. Silakan jalankan simulasi pada menu **Modul PCI** dan **Modul SDI** terlebih dahulu agar Dashboard Komparasi dapat ditampilkan.")

























