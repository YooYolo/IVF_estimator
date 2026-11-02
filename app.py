import os
import re
import zipfile
import xml.etree.ElementTree as ET
from itertools import combinations
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy import stats
import streamlit as st

# ---------------------------------------------------------------------------
# Page Configuration & Constants
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="IVF 统计与漏斗预测 / IVF Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

STUDY_KEY = "study"
COMMUNITY_KEY = "community"

DATASET_COLORS = {
    STUDY_KEY: "#000000",      # Black dots for Study dataset
    COMMUNITY_KEY: "#FF8C00",  # Orange dots for Community dataset
}

# Scatter marker styling: Study dots smaller with alpha=0.3; Community dots orange & distinct
DATASET_MARKER_STYLE = {
    STUDY_KEY: {"size": 5, "opacity": 0.3, "line_width": 0},
    COMMUNITY_KEY: {"size": 8, "opacity": 0.85, "line_width": 0.5},
}

STUDY_XLSX_PATH = os.path.join(
    os.path.dirname(__file__),
    "data",
    "data Diminished ovarian reserve is associated to euploidy rate.xlsx",
)

COMMUNITY_SHEET_ID = "1sHdyPNEOtVj8Yi_vtrNq4T2b3LaqdqYU"
COMMUNITY_GID = "920796724"
COMMUNITY_CSV_URL = (
    f"https://docs.google.com/spreadsheets/d/{COMMUNITY_SHEET_ID}"
    f"/export?format=csv&gid={COMMUNITY_GID}"
)

# ---------------------------------------------------------------------------
# Succinct Bilingual Labels (Chinese & English)
# ---------------------------------------------------------------------------
PARAM_LABELS_I18N: Dict[str, Dict[str, str]] = {
    "zh": {
        "age": "年龄",
        "amh": "AMH",
        "eggs_retrieved": "取卵数",
        "eggs_mature": "成熟卵数",
        "fertilized": "受精数",
        "blastocysts": "囊胚数",
        "euploid": "PGT正常数",
        "mature_rate": "成熟率",
        "fertilization_rate": "受精率",
        "blastocyst_rate": "养囊率",
        "euploid_rate": "过检率",
        "blast_per_egg_rate": "囊胚/取卵率",
        "euploid_per_egg_rate": "过检/取卵率",
        "gonadotropin_total_iu": "促排总剂量(IU)",
        "fsh_baseline": "基础FSH",
        "afc": "基础卵泡(AFC)",
        "bmi": "BMI",
        "estradiol": "雌二醇",
        "follicles_aspirated": "穿刺卵泡数",
        "aneuploid": "异常胚胎数",
        "fsh_dose_iu": "FSH剂量",
        "menopur_dose_iu": "美诺孕剂量",
        "mosaic": "嵌合体数",
        "retrieval_round": "促排周期",
    },
    "en": {
        "age": "Age",
        "amh": "AMH",
        "eggs_retrieved": "Eggs Retrieved",
        "eggs_mature": "Eggs Matured",
        "fertilized": "Fertilized",
        "blastocysts": "Blastocysts",
        "euploid": "Euploid",
        "mature_rate": "Mature Rate",
        "fertilization_rate": "Fertilize Rate",
        "blastocyst_rate": "Blastocyst Rate",
        "euploid_rate": "Euploid Rate",
        "blast_per_egg_rate": "Blast / Egg",
        "euploid_per_egg_rate": "Euploid / Egg",
        "gonadotropin_total_iu": "FSH+Menopur (IU)",
        "fsh_baseline": "Baseline FSH",
        "afc": "AFC",
        "bmi": "BMI",
        "estradiol": "Estradiol",
        "follicles_aspirated": "Follicles",
        "aneuploid": "Aneuploid",
        "fsh_dose_iu": "FSH Dose",
        "menopur_dose_iu": "Menopur Dose",
        "mosaic": "Mosaic",
        "retrieval_round": "Round",
    },
}

UI_TEXT: Dict[str, Dict[str, str]] = {
    "zh": {
        "title": "🧬 IVF 统计与养囊漏斗预测看板",
        "subtitle": (
            "**⚫ 临床研究 (黑点):** [Frontiers 2025 论文数据](https://www.frontiersin.org/journals/endocrinology/articles/10.3389/fendo.2024.1535776/full) &nbsp;|&nbsp; "
            "**🟠 群友数据 (橙点):** [在线共享表格](https://docs.google.com/spreadsheets/d/1sHdyPNEOtVj8Yi_vtrNq4T2b3LaqdqYU/edit?gid=920796724#gid=920796724)"
        ),
        "study_name": "临床研究 (Study)",
        "comm_name": "群友数据 (Community)",
        "both_name": "全部显示 (研究 + 群友)",
        "sidebar_header": "⚙️ 看板设置",
        "display_mode": "显示数据集",
        "status_header": "📡 数据状态",
        "cycles_loaded": "个周期",
        "comm_offline": "🟠 无法直接连接 Google Sheet。",
        "refresh_btn": "🔄 刷新在线表格",
        "upload_exp": "📂 上传本地表格 (.csv/.xlsx)",
        "upload_label": "选择文件",
        "legend_text": "**图例:** ⚫ 临床研究 (小黑点 α=0.3) | 🟠 群友数据 (橙点)",
        "tab_knn": "🎯 1. 漏斗预测 (k-NN)",
        "tab_corr": "📈 2. 相关性分析",
        "tab_data": "🗂️ 3. 汇总统计与原始数据",
        "corr_header": "交互式散点图与相关性分析 (Benjamini–Hochberg 多重检验校正)",
        "x_param": "X 轴参数",
        "y_param": "Y 轴参数",
        "stats_ds": "统计所用数据集 (默认: 临床研究)",
        "compare_both": "并列对比两组",
        "same_param_info": "请选择两个不同的参数。",
        "insufficient_data": "有效配对样本不足。",
        "metric_n": "样本量 N",
        "metric_r": "Pearson r",
        "metric_p": "原始 p 值",
        "metric_q": "BH 校正 q 值",
        "metric_rho": "Spearman ρ",
        "col_sig": "显著性",
        "sig_ref_title": "**显著性标记说明 (`*` 基于 BH 校正 q 值):** &nbsp; `*` $q < 0.05$ &nbsp;|&nbsp; `**` $q < 0.01$ &nbsp;|&nbsp; `***` $q < 0.001$ &nbsp;|&nbsp; `****` $q < 0.0001$ &nbsp;|&nbsp; `ns` $q \\ge 0.05$ (不显著)",
        "show_trend": "显示线性趋势线",
        "all_pairs_header": "📋 所有参数对相关性汇总表 (BH 校正)",
        "filter_param": "按参数筛选",
        "all_params": "全部参数",
        "only_sig": "仅显示显著相关 (BH q < 0.05)",
        "col_p1": "参数 1",
        "col_p2": "参数 2",
        "knn_header": "🎯 最近邻 (k-NN) 促排漏斗预测",
        "knn_desc": "可自由选择匹配参数（默认 **年龄** 与 **AMH**，也可加入 **取卵数、成熟卵数、受精数、囊胚数** 等，默认 $k=10$）。可在下方表格取消勾选 **`Include`** 剔除特定邻居并实时重算。",
        "knn_match_params": "匹配参数 (可自由增删，如：受精数、囊胚数等)",
        "select_at_least_one_param": "请至少选择 1 个匹配参数。",
        "input_age": "年龄",
        "input_amh": "AMH (ng/mL)",
        "input_k": "最近邻数量 (k)",
        "knn_pool": "匹配样本库",
        "pool_study": "仅临床研究",
        "pool_comm": "仅群友数据",
        "pool_blend": "混合样本库 (研究 + 群友)",
        "adv_filter_exp": "🔍 可选筛选条件 & 已知取卵数折算",
        "req_complete": "仅匹配五阶段数据完整的周期",
        "min_eggs": "最低取卵数过滤",
        "anchor_eggs": "按我已知的取卵数折算漏斗",
        "my_eggs": "我的取卵数",
        "no_neighbors": "未找到匹配的最近邻样本。",
        "step1_neighbors": "#### 1. 匹配到的最近邻样本 (取消勾选 `Include` 可剔除)",
        "active_count": "当前使用 **{active} / {total}** 个最近邻（已剔除 {removed} 个）。",
        "all_unchecked": "请至少保留勾选 1 个最近邻样本。",
        "step2_funnel": "#### 2. 养囊与过检漏斗预测结果",
        "of_retrieved": "% 占取卵数",
        "funnel_chart_title": "漏斗折线图 (均值、中位数及个体轨迹)",
        "mean_line": "均值 (Mean)",
        "median_line": "中位数 (Median)",
        "map_chart_title": "最近邻分布图 ({x_label} vs {y_label})",
        "map_x_axis": "分布图 X 轴",
        "map_y_axis": "分布图 Y 轴",
        "your_input": "您的输入",
        "selected_nb": "入选最近邻",
        "funnel_table_title": "##### 📋 各阶段漏斗统计明细",
        "summary_header": "#### 数据集汇总统计",
        "summary_ds_select": "选择数据集",
        "summary_stat_select": "选择统计量",
        "raw_header": "#### 标准化周期数据明细",
    },
    "en": {
        "title": "🧬 IVF Statistics & Funnel Estimator",
        "subtitle": (
            "**⚫ Study (Black):** [Frontiers in Endocrinology 2025](https://www.frontiersin.org/journals/endocrinology/articles/10.3389/fendo.2024.1535776/full) &nbsp;|&nbsp; "
            "**🟠 Community (Orange):** [Live Google Sheet](https://docs.google.com/spreadsheets/d/1sHdyPNEOtVj8Yi_vtrNq4T2b3LaqdqYU/edit?gid=920796724#gid=920796724)"
        ),
        "study_name": "Study",
        "comm_name": "Community",
        "both_name": "Both (Study + Community)",
        "sidebar_header": "⚙️ Settings",
        "display_mode": "Display Dataset",
        "status_header": "📡 Data Status",
        "cycles_loaded": "cycles",
        "comm_offline": "🟠 Could not reach live Google Sheet.",
        "refresh_btn": "🔄 Refresh Google Sheet",
        "upload_exp": "📂 Upload Community CSV/XLSX",
        "upload_label": "Upload file (.csv or .xlsx)",
        "legend_text": "**Legend:** ⚫ Study (size=5, α=0.3) | 🟠 Community",
        "tab_knn": "🎯 1. Funnel Estimator (k-NN)",
        "tab_corr": "📈 2. Correlations",
        "tab_data": "🗂️ 3. Summary & Data",
        "corr_header": "Scatterplot & Correlation Analysis (Benjamini–Hochberg Adjusted)",
        "x_param": "X Parameter",
        "y_param": "Y Parameter",
        "stats_ds": "Stats Dataset (Default: Study)",
        "compare_both": "Both Side-by-Side",
        "same_param_info": "Select two different parameters.",
        "insufficient_data": "Insufficient paired data.",
        "metric_n": "N",
        "metric_r": "Pearson r",
        "metric_p": "Raw p",
        "metric_q": "BH Adj. q",
        "metric_rho": "Spearman ρ",
        "col_sig": "Sig.",
        "sig_ref_title": "**Significance Reference (`*` based on BH Adj. q):** &nbsp; `*` $q < 0.05$ &nbsp;|&nbsp; `**` $q < 0.01$ &nbsp;|&nbsp; `***` $q < 0.001$ &nbsp;|&nbsp; `****` $q < 0.0001$ &nbsp;|&nbsp; `ns` $q \\ge 0.05$",
        "show_trend": "Show Trendline",
        "all_pairs_header": "📋 All-Pairs Correlation Table (BH Adjusted)",
        "filter_param": "Filter by Parameter",
        "all_params": "All",
        "only_sig": "Only Significant (BH q < 0.05)",
        "col_p1": "Param 1",
        "col_p2": "Param 2",
        "knn_header": "🎯 Nearest-Neighbor (k-NN) Attrition Funnel",
        "knn_desc": "Customize matching parameters (defaults to **Age** and **AMH**; can add **Eggs Retrieved, Matured, Fertilized, Blastocysts**, etc., default $k=10$). Uncheck **`Include`** in the table to remove specific neighbors and recalculate.",
        "knn_match_params": "Matching Parameters (add/remove e.g. Fertilized, Blastocysts, Eggs Retrieved)",
        "select_at_least_one_param": "Please select at least 1 matching parameter.",
        "input_age": "Age",
        "input_amh": "AMH (ng/mL)",
        "input_k": "Neighbors (k)",
        "knn_pool": "Neighbor Pool",
        "pool_study": "Study Only",
        "pool_comm": "Community Only",
        "pool_blend": "Blend (Both)",
        "adv_filter_exp": "🔍 Optional Filters & Starting Eggs Override",
        "req_complete": "Require complete 5-stage funnel data",
        "min_eggs": "Min Eggs Retrieved",
        "anchor_eggs": "Scale funnel to my known Eggs Retrieved",
        "my_eggs": "My Eggs Retrieved",
        "no_neighbors": "No matching neighbors found.",
        "step1_neighbors": "#### 1. Nearest Neighbors (Uncheck `Include` to Remove)",
        "active_count": "Using **{active} / {total}** neighbors ({removed} removed).",
        "all_unchecked": "Please keep at least 1 neighbor checked.",
        "step2_funnel": "#### 2. Estimated Attrition Funnel",
        "of_retrieved": "% of retrieved",
        "funnel_chart_title": "Attrition Funnel (Mean, Median & Neighbors)",
        "mean_line": "Mean",
        "median_line": "Median",
        "map_chart_title": "Neighbors in ({x_label} vs {y_label}) Space",
        "map_x_axis": "Map X-Axis",
        "map_y_axis": "Map Y-Axis",
        "your_input": "Your Input",
        "selected_nb": "Neighbors",
        "funnel_table_title": "##### 📋 Funnel Stage Summary",
        "summary_header": "#### Summary Statistics by Dataset",
        "summary_ds_select": "Select Dataset",
        "summary_stat_select": "Select Statistic",
        "raw_header": "#### Harmonized Records",
    },
}

FUNNEL_STAGE_COLS = [
    "eggs_retrieved",
    "eggs_mature",
    "fertilized",
    "blastocysts",
    "euploid",
]

KNN_PARAM_SPECS: Dict[str, Dict[str, object]] = {
    "age": {"min": 18.0, "max": 50.0, "default": 36.0, "step": 0.5, "fmt": "%.1f", "unit": ""},
    "amh": {"min": 0.01, "max": 25.0, "default": 1.50, "step": 0.1, "fmt": "%.2f", "unit": " (ng/mL)"},
    "eggs_retrieved": {"min": 0.0, "max": 60.0, "default": 10.0, "step": 1.0, "fmt": "%.0f", "unit": ""},
    "eggs_mature": {"min": 0.0, "max": 60.0, "default": 8.0, "step": 1.0, "fmt": "%.0f", "unit": ""},
    "fertilized": {"min": 0.0, "max": 50.0, "default": 6.0, "step": 1.0, "fmt": "%.0f", "unit": ""},
    "blastocysts": {"min": 0.0, "max": 30.0, "default": 3.0, "step": 1.0, "fmt": "%.0f", "unit": ""},
    "euploid": {"min": 0.0, "max": 20.0, "default": 1.0, "step": 1.0, "fmt": "%.0f", "unit": ""},
    "gonadotropin_total_iu": {"min": 0.0, "max": 10000.0, "default": 2400.0, "step": 75.0, "fmt": "%.0f", "unit": ""},
    "afc": {"min": 0.0, "max": 60.0, "default": 10.0, "step": 1.0, "fmt": "%.0f", "unit": ""},
    "fsh_baseline": {"min": 0.1, "max": 50.0, "default": 7.5, "step": 0.5, "fmt": "%.1f", "unit": " (mIU/mL)"},
    "bmi": {"min": 14.0, "max": 50.0, "default": 22.0, "step": 0.5, "fmt": "%.1f", "unit": ""},
    "mature_rate": {"min": 0.0, "max": 1.0, "default": 0.80, "step": 0.05, "fmt": "%.2f", "unit": " (0–1)"},
    "fertilization_rate": {"min": 0.0, "max": 1.0, "default": 0.75, "step": 0.05, "fmt": "%.2f", "unit": " (0–1)"},
    "blastocyst_rate": {"min": 0.0, "max": 1.0, "default": 0.50, "step": 0.05, "fmt": "%.2f", "unit": " (0–1)"},
    "euploid_rate": {"min": 0.0, "max": 1.0, "default": 0.45, "step": 0.05, "fmt": "%.2f", "unit": " (0–1)"},
}

# Exact (case-insensitive) translations for IVF Reasons (Italian study values + common EN/ZH terms)
IVF_REASON_EXACT_MAP: Dict[str, Tuple[str, str]] = {
    # Italian values from Frontiers 2025 Study Dataset (pma_indicazione_des) -> (Chinese, English)
    "disovulatoria": ("排卵障碍", "Ovulatory Disorder"),
    "endometriosi": ("子宫内膜异位症", "Endometriosis"),
    "fattore femminile": ("女方因素", "Female Factor"),
    "fattore maschile": ("男方因素", "Male Factor"),
    "fattore maschile grave": ("严重男方因素", "Severe Male Factor"),
    "fattore maschile lieve": ("轻度男方因素", "Mild Male Factor"),
    "fattore misto maschile e femminile": ("男女双方混合因素", "Combined Male & Female Factor"),
    "fattore multiplo femminile": ("多发女方因素", "Multiple Female Factors"),
    "fattore tubarico parziale": ("部分输卵管因素", "Partial Tubal Factor"),
    "fattore tubarico/aderenziale": ("输卵管/盆腔粘连", "Tubal / Pelvic Adhesions"),
    "genetico": ("遗传因素", "Genetic Factor"),
    "idiopatica": ("不明原因", "Unexplained (Idiopathic)"),
    "poliabortività": ("复发性流产", "Recurrent Pregnancy Loss"),
    "ridotta riserva ovarica": ("卵巢储备功能减退 (DOR)", "Diminished Ovarian Reserve (DOR)"),
    # Common English terms -> (Chinese, English)
    "dor": ("卵巢储备功能减退 (DOR)", "DOR"),
    "diminished ovarian reserve": ("卵巢储备功能减退 (DOR)", "Diminished Ovarian Reserve"),
    "male factor": ("男方因素", "Male Factor"),
    "mfi": ("男方因素 (MFI)", "Male Factor (MFI)"),
    "female factor": ("女方因素", "Female Factor"),
    "tubal": ("输卵管因素", "Tubal Factor"),
    "tubal factor": ("输卵管因素", "Tubal Factor"),
    "pcos": ("多囊卵巢综合征 (PCOS)", "PCOS"),
    "endometriosis": ("子宫内膜异位症", "Endometriosis"),
    "unexplained": ("不明原因", "Unexplained"),
    "age": ("高龄因素", "Age Factor"),
    "ama": ("高龄因素 (AMA)", "Advanced Maternal Age"),
    "advanced maternal age": ("高龄因素", "Advanced Maternal Age"),
    "rpl": ("复发性流产 (RPL)", "Recurrent Pregnancy Loss"),
    "recurrent pregnancy loss": ("复发性流产", "Recurrent Pregnancy Loss"),
    "genetic": ("遗传/基因因素", "Genetic Factor"),
    "pgt-m": ("单基因遗传病 (PGT-M)", "PGT-M (Monogenic)"),
    "pgt-sr": ("染色体结构异常 (PGT-SR)", "PGT-SR (Structural Rearrangement)"),
    "egg freezing": ("冻卵/生育力保存", "Egg Freezing"),
    "fertility preservation": ("生育力保存", "Fertility Preservation"),
}

# Substring/phrase mapping for free-text Chinese <-> English community reasons
ZH_TO_EN_REASON_TERMS: List[Tuple[str, str]] = [
    ("卵巢储备功能减退", "DOR"),
    ("卵巢功能减退", "DOR"),
    ("卵巢早衰", "POI/DOR"),
    ("早衰", "DOR"),
    ("多囊卵巢", "PCOS"),
    ("多囊", "PCOS"),
    ("子宫内膜异位症", "Endometriosis"),
    ("子宫内膜异位", "Endometriosis"),
    ("内异症", "Endometriosis"),
    ("巧囊", "Endometrioma"),
    ("腺肌症", "Adenomyosis"),
    ("双侧输卵管切除", "Bilateral Salpingectomy"),
    ("输卵管堵塞", "Tubal Blockage"),
    ("输卵管积水", "Hydrosalpinx"),
    ("输卵管粘连", "Tubal Adhesions"),
    ("输卵管", "Tubal Factor"),
    ("宫外孕", "Ectopic History"),
    ("排卵障碍", "Ovulatory Disorder"),
    ("少弱畸精", "OAT (Male Factor)"),
    ("少弱精", "Oligoasthenozoospermia"),
    ("弱精", "Asthenozoospermia (Male)"),
    ("少精", "Oligozoospermia (Male)"),
    ("畸精", "Teratozoospermia (Male)"),
    ("无精", "Azoospermia (Male)"),
    ("精子碎片率高", "High Sperm DFI"),
    ("碎片率高", "High DFI"),
    ("精索静脉曲张", "Varicocele"),
    ("男方因素", "Male Factor"),
    ("男方", "Male Factor"),
    ("女方因素", "Female Factor"),
    ("复发性流产", "Recurrent Pregnancy Loss"),
    ("习惯性流产", "Recurrent Miscarriage"),
    ("胎停", "Missed Miscarriage"),
    ("生化", "Biochemical Pregnancy"),
    ("流产", "Miscarriage"),
    ("平衡易位", "Balanced Translocation"),
    ("罗氏易位", "Robertsonian Translocation"),
    ("染色体异常", "Chromosomal Abnormality"),
    ("染色体", "Chromosomal"),
    ("单基因病", "Monogenic (PGT-M)"),
    ("单基因", "Monogenic (PGT-M)"),
    ("遗传", "Genetic"),
    ("基因", "Genetic"),
    ("不明原因不孕", "Unexplained Infertility"),
    ("不明原因", "Unexplained"),
    ("高龄", "Age Factor"),
    ("冻卵", "Egg Freezing"),
    ("选性别", "Family Balancing"),
    ("二胎", "Second Child"),
    ("三代", "PGT"),
]

PROTOCOL_EXACT_MAP: Dict[str, Tuple[str, str]] = {
    "antagonista fisso": ("固定拮抗剂方案", "Fixed Antagonist"),
    "antagonista variabile": ("灵活拮抗剂方案", "Flexible Antagonist"),
    "flare": ("微刺激/短效激发 (Flare)", "Flare-up"),
    "long": ("长方案", "Long Agonist"),
    "short": ("短方案", "Short Agonist"),
    "spontaneo": ("自然周期", "Natural Cycle"),
}


def translate_ivf_reason(val: object, lang: str) -> str:
    if val is None or pd.isna(val):
        return ""
    s = str(val).strip()
    if not s or s.lower() == "nan":
        return ""

    key = s.lower()
    if key in IVF_REASON_EXACT_MAP:
        zh_val, en_val = IVF_REASON_EXACT_MAP[key]
        return zh_val if lang == "zh" else en_val

    if lang == "en" and re.search(r"[\u4e00-\u9fff]", s):
        out = s
        for zh_term, en_term in ZH_TO_EN_REASON_TERMS:
            if zh_term in out:
                out = out.replace(zh_term, f" {en_term} ")
        out = (
            out.replace("，", ", ")
            .replace("、", ", ")
            .replace("；", "; ")
            .replace("＋", " + ")
        )
        out = re.sub(r"\s+", " ", out).strip(" ,;")
        return out

    return s


def translate_protocol(val: object, lang: str) -> str:
    if val is None or pd.isna(val):
        return ""
    s = str(val).strip()
    if not s or s.lower() == "nan":
        return ""
    key = s.lower()
    if key in PROTOCOL_EXACT_MAP:
        zh_val, en_val = PROTOCOL_EXACT_MAP[key]
        return zh_val if lang == "zh" else en_val
    return s


# ---------------------------------------------------------------------------
# Helper Functions for Parsing & Cleaning
# ---------------------------------------------------------------------------
def _parse_numeric_series(series: pd.Series, is_rate: bool = False) -> pd.Series:
    if series is None:
        return pd.Series(dtype=float)
    if pd.api.types.is_numeric_dtype(series):
        out = pd.to_numeric(series, errors="coerce").astype(float)
    else:
        s = series.astype(str).str.strip()
        has_pct = s.str.contains("%", na=False)
        cleaned = (
            s.str.replace("%", "", regex=False)
            .str.replace(",", ".", regex=False)
            .str.extract(r"([-+]?\d*\.?\d+)", expand=False)
        )
        out = pd.to_numeric(cleaned, errors="coerce").astype(float)
        if is_rate:
            out = np.where(has_pct | (out > 1.5), out / 100.0, out)
            out = pd.Series(out, index=series.index, dtype=float)
    if is_rate:
        out = out.mask((out < 0) | (out > 1.5))
        out = out.clip(lower=0.0, upper=1.0)
    return out


def _safe_ratio(num: pd.Series, den: pd.Series) -> pd.Series:
    ratio = np.where((den > 0) & num.notna() & den.notna(), num / den, np.nan)
    s = pd.Series(ratio, index=num.index, dtype=float)
    return s.mask((s < 0) | (s > 1.2)).clip(0.0, 1.0)


def _read_xlsx_fallback(xlsx_path: str) -> pd.DataFrame:
    def col_to_idx(col_str: str) -> int:
        idx = 0
        for ch in col_str:
            idx = idx * 26 + (ord(ch) - ord("A") + 1)
        return idx - 1

    with zipfile.ZipFile(xlsx_path) as z:
        ss = []
        if "xl/sharedStrings.xml" in z.namelist():
            ss_xml = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in ss_xml.iter():
                if si.tag.endswith("si"):
                    ss.append("".join(t.text or "" for t in si.iter() if t.tag.endswith("t")))
        sheet1 = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
        rows = []
        for row in sheet1.iter():
            if row.tag.endswith("row"):
                r_dict = {}
                for c in row:
                    ref = c.attrib.get("r", "")
                    m = re.match(r"[A-Z]+", ref)
                    if not m:
                        continue
                    idx = col_to_idx(m.group(0))
                    t = c.attrib.get("t")
                    v = next((ch.text for ch in c if ch.tag.endswith("v")), None)
                    if t == "s" and v is not None:
                        v = ss[int(v)]
                    r_dict[idx] = v
                rows.append(r_dict)
        max_col = max(max(r.keys()) for r in rows if r)
        headers = [rows[0].get(i, f"Unnamed: {i}") for i in range(max_col + 1)]
        data = [[r.get(i, None) for i in range(max_col + 1)] for r in rows[1:]]
        return pd.DataFrame(data, columns=headers)


# ---------------------------------------------------------------------------
# Dataset Loaders
# ---------------------------------------------------------------------------
@st.cache_data(show_spinner=False)
def load_study_dataset(xlsx_path: str = STUDY_XLSX_PATH) -> pd.DataFrame:
    if not os.path.exists(xlsx_path):
        return pd.DataFrame()

    try:
        raw = pd.read_excel(xlsx_path, engine="openpyxl")
    except Exception:
        raw = _read_xlsx_fallback(xlsx_path)

    first_col = raw.columns[0]
    mask_valid_id = raw[first_col].astype(str).str.strip().str.startswith("ICH", na=False)
    raw = raw.loc[mask_valid_id].copy().reset_index(drop=True)

    df = pd.DataFrame()
    df["record_id"] = raw[first_col].astype(str).str.strip()
    df["dataset_key"] = STUDY_KEY

    df["age"] = _parse_numeric_series(raw.get("paz_eta_donna"))
    df["amh"] = _parse_numeric_series(raw.get("esa_last_amh"))
    df["fsh_baseline"] = _parse_numeric_series(raw.get("esa_last_fsh"))
    df["afc"] = _parse_numeric_series(raw.get("ecg_afc"))
    df["bmi"] = _parse_numeric_series(raw.get("ana_bmi_donna"))
    df["estradiol"] = _parse_numeric_series(raw.get("estradiolo"))
    df["gonadotropin_total_iu"] = _parse_numeric_series(raw.get("mon_gonad_dosaggio"))

    df["follicles_aspirated"] = _parse_numeric_series(raw.get("foll_aspirati"))
    df["eggs_retrieved"] = _parse_numeric_series(raw.get("ovo_prelevati"))
    df["eggs_mature"] = _parse_numeric_series(raw.get("ovo_post_m2"))

    fert_2pn = _parse_numeric_series(raw.get("ovo_2pn"))
    fert_tot = _parse_numeric_series(raw.get("ovo_fecondati"))
    df["fertilized"] = fert_2pn.fillna(fert_tot)

    df["blastocysts"] = _parse_numeric_series(raw.get("emb_biopsati"))
    df["euploid"] = _parse_numeric_series(raw.get("emb_biopsia_pgs_euploide"))
    df["aneuploid"] = _parse_numeric_series(raw.get("emb_biopsia_pgs_aneuploide"))

    df["fsh_dose_iu"] = np.nan
    df["menopur_dose_iu"] = np.nan
    df["mosaic"] = np.nan
    df["retrieval_round"] = np.nan

    df["mature_rate"] = _safe_ratio(df["eggs_mature"], df["eggs_retrieved"])
    df["fertilization_rate"] = _safe_ratio(df["fertilized"], df["eggs_mature"])
    df["blastocyst_rate"] = _safe_ratio(df["blastocysts"], df["fertilized"])
    df["euploid_rate"] = _safe_ratio(df["euploid"], df["blastocysts"])
    df["blast_per_egg_rate"] = _safe_ratio(df["blastocysts"], df["eggs_retrieved"])
    df["euploid_per_egg_rate"] = _safe_ratio(df["euploid"], df["eggs_retrieved"])

    df["ivf_reason"] = raw.get("pma_indicazione_des", pd.Series("", index=raw.index)).fillna("").astype(str).str.strip()
    df["protocol"] = raw.get("mon_protocollo_des", pd.Series("", index=raw.index)).fillna("").astype(str).str.strip()
    df["trigger"] = raw.get("mon_hcg_des", pd.Series("", index=raw.index)).fillna("").astype(str).str.strip()
    df["note"] = ""

    return df


def _find_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    norm_map = {re.sub(r"\s+", " ", str(c).strip().lower()): c for c in df.columns}
    for cand in candidates:
        key = re.sub(r"\s+", " ", cand.strip().lower())
        if key in norm_map:
            return norm_map[key]
    return None


def harmonize_community_df(raw: pd.DataFrame) -> pd.DataFrame:
    if raw is None or raw.empty:
        return pd.DataFrame()

    c_name = _find_col(raw, ["Name", "ID", "Patient"])
    c_round = _find_col(raw, ["Egg retrieval round", "Retrieval round", "Round"])
    c_age = _find_col(raw, ["age", "Age"])
    c_amh = _find_col(raw, ["AMH", "amh"])
    c_fsh = _find_col(raw, ["FSH", "fsh"])
    c_meno = _find_col(raw, ["Menopur", "menopur"])
    c_fsh_meno = _find_col(raw, ["FSH+Menopur", "FSH + Menopur"])
    c_ret = _find_col(raw, ["Egg retreived", "Egg retrieved", "Eggs retrieved"])
    c_mat = _find_col(raw, ["Egg matured", "Eggs matured", "Mature", "MII"])
    c_fert = _find_col(raw, ["fertilized", "Fertilized", "2PN"])
    c_blast = _find_col(raw, ["Blastocyst", "Blastocysts"])
    c_eu = _find_col(raw, ["Euploid", "Euploids"])
    c_mat_rate = _find_col(raw, ["Mature rate", "Maturation rate"])
    c_fert_rate = _find_col(raw, ["Fertilize rate", "Fertilization rate"])
    c_blast_rate = _find_col(raw, ["Blastocyst rate", "Blast rate"])
    c_eu_rate = _find_col(raw, ["Eu rate", "Euploid rate", "Euploidy rate"])
    c_blast_egg = _find_col(raw, ["Blastocyst/egg", "Blast/egg"])
    c_mosaic = _find_col(raw, ["Mosaic", "Mosaics"])
    c_trig = _find_col(raw, ["Trigger"])
    c_reason = _find_col(raw, ["IVF reason", "Reason", "Diagnosis"])
    c_note = _find_col(raw, ["Note", "Notes"])

    df = pd.DataFrame(index=raw.index)
    name_s = raw[c_name].astype(str).str.strip() if c_name else pd.Series(raw.index.astype(str), index=raw.index)
    round_num = _parse_numeric_series(raw[c_round]) if c_round else pd.Series(np.nan, index=raw.index)

    rec_ids = []
    for idx_val, nm, rn in zip(raw.index, name_s, round_num):
        nm_clean = nm if (nm and nm.lower() != "nan") else f"Row-{idx_val + 1}"
        if nm_clean.endswith(".0"):
            nm_clean = nm_clean[:-2]
        if pd.notna(rn):
            rec_ids.append(f"#{nm_clean}-R{int(rn)}")
        else:
            rec_ids.append(f"#{nm_clean}")

    df["record_id"] = rec_ids
    df["dataset_key"] = COMMUNITY_KEY
    df["retrieval_round"] = round_num
    df["age"] = _parse_numeric_series(raw[c_age]) if c_age else np.nan
    df["amh"] = _parse_numeric_series(raw[c_amh]) if c_amh else np.nan

    df["fsh_dose_iu"] = _parse_numeric_series(raw[c_fsh]) if c_fsh else np.nan
    df["menopur_dose_iu"] = _parse_numeric_series(raw[c_meno]) if c_meno else np.nan
    fsh_plus_meno = _parse_numeric_series(raw[c_fsh_meno]) if c_fsh_meno else pd.Series(np.nan, index=raw.index)
    summed_gonad = df["fsh_dose_iu"].fillna(0) + df["menopur_dose_iu"].fillna(0)
    summed_gonad = summed_gonad.where(df["fsh_dose_iu"].notna() | df["menopur_dose_iu"].notna(), np.nan)
    df["gonadotropin_total_iu"] = fsh_plus_meno.fillna(summed_gonad)

    df["eggs_retrieved"] = _parse_numeric_series(raw[c_ret]) if c_ret else np.nan
    df["eggs_mature"] = _parse_numeric_series(raw[c_mat]) if c_mat else np.nan
    df["fertilized"] = _parse_numeric_series(raw[c_fert]) if c_fert else np.nan
    df["blastocysts"] = _parse_numeric_series(raw[c_blast]) if c_blast else np.nan
    df["euploid"] = _parse_numeric_series(raw[c_eu]) if c_eu else np.nan
    df["mosaic"] = _parse_numeric_series(raw[c_mosaic]) if c_mosaic else np.nan

    df["fsh_baseline"] = np.nan
    df["afc"] = np.nan
    df["bmi"] = np.nan
    df["estradiol"] = np.nan
    df["follicles_aspirated"] = np.nan
    df["aneuploid"] = np.nan

    calc_mat_rate = _safe_ratio(df["eggs_mature"], df["eggs_retrieved"])
    sheet_mat_rate = _parse_numeric_series(raw[c_mat_rate], is_rate=True) if c_mat_rate else pd.Series(np.nan, index=raw.index)
    df["mature_rate"] = calc_mat_rate.fillna(sheet_mat_rate)

    calc_fert_rate = _safe_ratio(df["fertilized"], df["eggs_mature"])
    sheet_fert_rate = _parse_numeric_series(raw[c_fert_rate], is_rate=True) if c_fert_rate else pd.Series(np.nan, index=raw.index)
    df["fertilization_rate"] = calc_fert_rate.fillna(sheet_fert_rate)

    calc_blast_rate = _safe_ratio(df["blastocysts"], df["fertilized"])
    sheet_blast_rate = _parse_numeric_series(raw[c_blast_rate], is_rate=True) if c_blast_rate else pd.Series(np.nan, index=raw.index)
    df["blastocyst_rate"] = calc_blast_rate.fillna(sheet_blast_rate)

    calc_eu_rate = _safe_ratio(df["euploid"], df["blastocysts"])
    sheet_eu_rate = _parse_numeric_series(raw[c_eu_rate], is_rate=True) if c_eu_rate else pd.Series(np.nan, index=raw.index)
    df["euploid_rate"] = calc_eu_rate.fillna(sheet_eu_rate)

    calc_blast_egg = _safe_ratio(df["blastocysts"], df["eggs_retrieved"])
    sheet_blast_egg = _parse_numeric_series(raw[c_blast_egg], is_rate=True) if c_blast_egg else pd.Series(np.nan, index=raw.index)
    df["blast_per_egg_rate"] = calc_blast_egg.fillna(sheet_blast_egg)

    df["euploid_per_egg_rate"] = _safe_ratio(df["euploid"], df["eggs_retrieved"])

    df["ivf_reason"] = raw[c_reason].fillna("").astype(str).str.strip() if c_reason else ""
    df["protocol"] = ""
    df["trigger"] = raw[c_trig].fillna("").astype(str).str.strip() if c_trig else ""
    df["note"] = raw[c_note].fillna("").astype(str).str.strip() if c_note else ""

    valid_mask = df[["age", "amh", "eggs_retrieved", "blastocysts", "euploid"]].notna().any(axis=1)
    valid_mask &= df["age"].isna() | ((df["age"] >= 18) & (df["age"] <= 55))
    return df.loc[valid_mask].reset_index(drop=True)


@st.cache_data(ttl=600, show_spinner=False)
def load_community_dataset_live(csv_url: str = COMMUNITY_CSV_URL) -> Tuple[pd.DataFrame, Optional[str]]:
    try:
        raw = pd.read_csv(csv_url)
        df = harmonize_community_df(raw)
        return df, None
    except Exception as exc:
        return pd.DataFrame(), str(exc)


# ---------------------------------------------------------------------------
# Statistical & Correlation Analysis with Benjamini-Hochberg Adjustment
# ---------------------------------------------------------------------------
def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    p_arr = np.asarray(p_values, dtype=float)
    q_arr = np.full_like(p_arr, np.nan, dtype=float)
    valid_mask = ~np.isnan(p_arr)
    m = int(valid_mask.sum())
    if m == 0:
        return q_arr

    p_valid = p_arr[valid_mask]
    order = np.argsort(p_valid)
    p_sorted = p_valid[order]

    ranks = np.arange(1, m + 1, dtype=float)
    q_sorted = p_sorted * (m / ranks)
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0.0, 1.0)

    q_valid = np.empty_like(q_sorted)
    q_valid[order] = q_sorted
    q_arr[valid_mask] = q_valid
    return q_arr


def sig_stars(p_val: float) -> str:
    """Return standard significance stars for a given p/q value."""
    if pd.isna(p_val):
        return ""
    if p_val < 0.0001:
        return "****"
    if p_val < 0.001:
        return "***"
    if p_val < 0.01:
        return "**"
    if p_val < 0.05:
        return "*"
    return "ns"


@st.cache_data(show_spinner=False)
def compute_all_pairs_correlations(df: pd.DataFrame, min_n: int = 8) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame()

    candidate_cols = [c for c in PARAM_LABELS_I18N["en"].keys() if c in df.columns]
    valid_cols = []
    for c in candidate_cols:
        s = df[c].dropna()
        if len(s) >= min_n and s.nunique() > 1:
            valid_cols.append(c)

    records = []
    for col_x, col_y in combinations(valid_cols, 2):
        sub = df[[col_x, col_y]].dropna()
        n = len(sub)
        if n < min_n or sub[col_x].nunique() < 2 or sub[col_y].nunique() < 2:
            continue

        x_vals = sub[col_x].to_numpy(dtype=float)
        y_vals = sub[col_y].to_numpy(dtype=float)

        try:
            r_pearson, p_pearson = stats.pearsonr(x_vals, y_vals)
        except Exception:
            r_pearson, p_pearson = np.nan, np.nan

        try:
            r_spearman, p_spearman = stats.spearmanr(x_vals, y_vals)
        except Exception:
            r_spearman, p_spearman = np.nan, np.nan

        records.append(
            {
                "param_x": col_x,
                "param_y": col_y,
                "n": n,
                "pearson_r": float(r_pearson),
                "pearson_p": float(p_pearson),
                "spearman_rho": float(r_spearman),
                "spearman_p": float(p_spearman),
            }
        )

    if not records:
        return pd.DataFrame()

    res = pd.DataFrame(records)
    res["pearson_q_bh"] = benjamini_hochberg(res["pearson_p"].to_numpy())
    res["spearman_q_bh"] = benjamini_hochberg(res["spearman_p"].to_numpy())
    res["abs_pearson_r"] = res["pearson_r"].abs()
    res = res.sort_values(by=["pearson_q_bh", "abs_pearson_r"], ascending=[True, False]).reset_index(drop=True)
    return res


def lookup_pair_stats(corr_df: pd.DataFrame, param_x: str, param_y: str) -> Optional[Dict]:
    if corr_df is None or corr_df.empty or param_x == param_y:
        return None
    m1 = (corr_df["param_x"] == param_x) & (corr_df["param_y"] == param_y)
    if m1.any():
        return corr_df.loc[m1].iloc[0].to_dict()
    m2 = (corr_df["param_x"] == param_y) & (corr_df["param_y"] == param_x)
    if m2.any():
        row = corr_df.loc[m2].iloc[0].to_dict()
        row["param_x"], row["param_y"] = param_x, param_y
        return row
    return None


# ---------------------------------------------------------------------------
# Nearest-Neighbor (k-NN) Attrition Funnel Functions
# ---------------------------------------------------------------------------
def find_nearest_neighbors(
    pool_df: pd.DataFrame,
    user_params: Optional[Dict[str, float]] = None,
    k: int = 10,
    user_age: Optional[float] = None,
    user_amh: Optional[float] = None,
) -> pd.DataFrame:
    if pool_df is None or pool_df.empty:
        return pd.DataFrame()

    target_params: Dict[str, float] = {}
    if user_params:
        for col, val in user_params.items():
            if col in pool_df.columns and val is not None and pd.notna(val):
                target_params[col] = float(val)
    else:
        if user_age is not None and "age" in pool_df.columns:
            target_params["age"] = float(user_age)
        if user_amh is not None and "amh" in pool_df.columns:
            target_params["amh"] = float(user_amh)

    feature_cols = list(target_params.keys())
    if not feature_cols:
        return pd.DataFrame()

    valid = pool_df.dropna(subset=feature_cols).copy()
    if valid.empty:
        return pd.DataFrame()

    sq_dist = np.zeros(len(valid), dtype=float)
    diff_cols: List[str] = []

    for col in feature_cols:
        u_val = target_params[col]
        c_mean = valid[col].mean()
        c_std = valid[col].std(ddof=0)
        c_std = c_std if (pd.notna(c_std) and c_std > 0) else 1.0

        z_col = (valid[col] - c_mean) / c_std
        u_z = (u_val - c_mean) / c_std
        sq_dist += (z_col - u_z) ** 2

        diff_col = f"_diff_{col}"
        valid[diff_col] = (valid[col] - u_val).abs()
        diff_cols.append(diff_col)

    valid["distance"] = np.sqrt(sq_dist)
    if "age" in target_params:
        valid["age_diff"] = valid["_diff_age"]
    if "amh" in target_params:
        valid["amh_diff"] = valid["_diff_amh"]

    valid = valid.sort_values(by=["distance"] + diff_cols, ascending=True).reset_index(drop=True)
    valid = valid.drop(columns=diff_cols, errors="ignore")
    neighbors = valid.head(int(k)).copy()
    neighbors.insert(0, "rank", np.arange(1, len(neighbors) + 1))
    return neighbors


def summarize_funnel(
    neighbors_df: pd.DataFrame,
    lang: str,
    starting_eggs_override: Optional[float] = None,
) -> pd.DataFrame:
    p_labels = PARAM_LABELS_I18N[lang]
    rows = []
    prev_mean = None
    prev_median = None

    for stage_col in FUNNEL_STAGE_COLS:
        stage_label = p_labels[stage_col]
        s = neighbors_df[stage_col].dropna() if stage_col in neighbors_df.columns else pd.Series(dtype=float)
        n_valid = len(s)
        mean_val = float(s.mean()) if n_valid > 0 else np.nan
        med_val = float(s.median()) if n_valid > 0 else np.nan
        q25_val = float(s.quantile(0.25)) if n_valid > 0 else np.nan
        q75_val = float(s.quantile(0.75)) if n_valid > 0 else np.nan
        min_val = float(s.min()) if n_valid > 0 else np.nan
        max_val = float(s.max()) if n_valid > 0 else np.nan

        step_retention_mean = (
            (mean_val / prev_mean)
            if (prev_mean and prev_mean > 0 and pd.notna(mean_val))
            else (1.0 if stage_col == "eggs_retrieved" else np.nan)
        )
        step_retention_med = (
            (med_val / prev_median)
            if (prev_median and prev_median > 0 and pd.notna(med_val))
            else (1.0 if stage_col == "eggs_retrieved" else np.nan)
        )

        prev_mean = mean_val
        prev_median = med_val

        rows.append(
            {
                "stage_col": stage_col,
                "Stage": stage_label,
                "N": n_valid,
                "Mean": mean_val,
                "Median": med_val,
                "Q25": q25_val,
                "Q75": q75_val,
                "Min": min_val,
                "Max": max_val,
                "Step % (Mean)": step_retention_mean * 100.0 if pd.notna(step_retention_mean) else np.nan,
                "Step % (Median)": step_retention_med * 100.0 if pd.notna(step_retention_med) else np.nan,
            }
        )

    summary = pd.DataFrame(rows)
    base_mean = summary.loc[summary["stage_col"] == "eggs_retrieved", "Mean"].iloc[0]
    if pd.notna(base_mean) and base_mean > 0:
        summary["Cumul. %"] = (summary["Mean"] / base_mean) * 100.0
    else:
        summary["Cumul. %"] = np.nan

    if starting_eggs_override is not None and starting_eggs_override > 0 and pd.notna(base_mean) and base_mean > 0:
        summary["Scaled Est."] = (summary["Mean"] / base_mean) * float(starting_eggs_override)

    return summary


# ---------------------------------------------------------------------------
# Wide Summary Table Builder
# ---------------------------------------------------------------------------
def build_wide_summary_table(
    df: pd.DataFrame,
    lang: str,
    ds_names: Dict[str, str],
    selected_ds_keys: List[str],
    selected_stats: List[str],
) -> pd.DataFrame:
    """
    Build a wide-format summary statistics table filtered by selected dataset(s)
    and selected statistic(s) (mean, median, std).
    """
    summary_cols = [
        "age",
        "amh",
        "gonadotropin_total_iu",
        "eggs_retrieved",
        "eggs_mature",
        "fertilized",
        "blastocysts",
        "euploid",
        "mature_rate",
        "fertilization_rate",
        "blastocyst_rate",
        "euploid_rate",
    ]
    present_cols = [c for c in summary_cols if c in df.columns]
    p_labels = PARAM_LABELS_I18N[lang]

    stat_names = (
        {"mean": "均值 (Mean)", "median": "中位数 (Median)", "std": "标准差 (SD)"}
        if lang == "zh"
        else {"mean": "Mean", "median": "Median", "std": "SD"}
    )

    records = []
    for ds_key in selected_ds_keys:
        sub = df[df["dataset_key"] == ds_key]
        if sub.empty:
            continue
        ds_label = ds_names[ds_key]
        agg_df = sub[present_cols].agg(["mean", "median", "std"])
        for stat_k in selected_stats:
            if stat_k not in agg_df.index:
                continue
            row_dict = {
                ("数据集" if lang == "zh" else "Dataset"): ds_label,
                ("统计量" if lang == "zh" else "Stat"): stat_names[stat_k],
            }
            for col in present_cols:
                row_dict[p_labels[col]] = agg_df.loc[stat_k, col]
            records.append(row_dict)

    wide_df = pd.DataFrame(records)
    if not wide_df.empty:
        index_cols = ["数据集", "统计量"] if lang == "zh" else ["Dataset", "Stat"]
        wide_df = wide_df.set_index(index_cols)
    return wide_df


# ---------------------------------------------------------------------------
# Main Dashboard UI
# ---------------------------------------------------------------------------
def main() -> None:
    study_df = load_study_dataset()
    comm_df, comm_err = load_community_dataset_live()

    # Sidebar Controls (Default Language: Chinese "zh")
    with st.sidebar:
        lang_choice = st.radio(
            "🌐 语言 / Language",
            options=["zh", "en"],
            format_func=lambda x: "中文" if x == "zh" else "English",
            index=0,
            horizontal=True,
        )
        t = UI_TEXT[lang_choice]
        p_labels = PARAM_LABELS_I18N[lang_choice]
        ds_names = {
            STUDY_KEY: t["study_name"],
            COMMUNITY_KEY: t["comm_name"],
        }

        st.header(t["sidebar_header"])

        display_mode = st.radio(
            t["display_mode"],
            options=["both", STUDY_KEY, COMMUNITY_KEY],
            format_func=lambda m: (
                t["both_name"]
                if m == "both"
                else (t["study_name"] if m == STUDY_KEY else t["comm_name"])
            ),
            index=0,
        )

        st.markdown("---")
        st.subheader(t["status_header"])
        st.caption(f"⚫ **{t['study_name']}:** {len(study_df):,} {t['cycles_loaded']}")

        if not comm_df.empty:
            st.caption(f"🟠 **{t['comm_name']}:** {len(comm_df):,} {t['cycles_loaded']}")
        else:
            st.warning(t["comm_offline"])
            if comm_err:
                st.caption(f"`{comm_err[:100]}`")

        if st.button(t["refresh_btn"]):
            load_community_dataset_live.clear()
            st.rerun()

        with st.expander(t["upload_exp"]):
            uploaded_comm = st.file_uploader(t["upload_label"], type=["csv", "xlsx"])
            if uploaded_comm is not None:
                try:
                    raw_up = (
                        pd.read_csv(uploaded_comm)
                        if uploaded_comm.name.endswith(".csv")
                        else pd.read_excel(uploaded_comm)
                    )
                    comm_df = harmonize_community_df(raw_up)
                    st.success(f"Loaded {len(comm_df):,} rows!")
                except Exception as up_e:
                    st.error(f"Error: {up_e}")

        st.markdown("---")
        st.caption(t["legend_text"])

    # Main Title
    st.title(t["title"])
    st.markdown(t["subtitle"])

    # Filter active dataframe
    frames = []
    if display_mode in ("both", STUDY_KEY) and not study_df.empty:
        frames.append(study_df)
    if display_mode in ("both", COMMUNITY_KEY) and not comm_df.empty:
        frames.append(comm_df)
    active_df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    # Precompute correlations with BH adjustment separately for each dataset
    study_corr_df = compute_all_pairs_correlations(study_df)
    comm_corr_df = compute_all_pairs_correlations(comm_df)

    # Default first page: Funnel Estimator (tab_knn)
    tab_knn, tab_corr, tab_data = st.tabs([t["tab_knn"], t["tab_corr"], t["tab_data"]])

    # =========================================================================
    # TAB 1 (Default): Nearest-Neighbor (k-NN) Attrition Funnel Estimator
    # =========================================================================
    with tab_knn:
        st.subheader(t["knn_header"])
        st.markdown(t["knn_desc"])

        top_c1, top_c2, top_c3 = st.columns([2.2, 1.0, 1.4])
        with top_c1:
            knn_param_options = list(KNN_PARAM_SPECS.keys())
            selected_knn_params = st.multiselect(
                t["knn_match_params"],
                options=knn_param_options,
                default=["age", "amh"],
                format_func=lambda k: p_labels.get(k, k),
            )
        with top_c2:
            k_neighbors = st.slider(t["input_k"], min_value=3, max_value=50, value=10, step=1)
        with top_c3:
            knn_pool_choice = st.radio(
                t["knn_pool"],
                options=[STUDY_KEY, COMMUNITY_KEY, "blend"],
                format_func=lambda x: (
                    t["pool_study"]
                    if x == STUDY_KEY
                    else (t["pool_comm"] if x == COMMUNITY_KEY else t["pool_blend"])
                ),
                index=2,
                horizontal=True,
            )

        user_params: Dict[str, float] = {}
        if selected_knn_params:
            n_cols = min(len(selected_knn_params), 4)
            param_cols = st.columns(n_cols)
            for idx, p_col in enumerate(selected_knn_params):
                spec = KNN_PARAM_SPECS[p_col]
                unit_str = str(spec.get("unit", ""))
                inp_label = f"{p_labels.get(p_col, p_col)}{unit_str}"
                with param_cols[idx % n_cols]:
                    user_params[p_col] = float(
                        st.number_input(
                            inp_label,
                            min_value=float(spec["min"]),
                            max_value=float(spec["max"]),
                            value=float(spec["default"]),
                            step=float(spec["step"]),
                            format=str(spec["fmt"]),
                            key=f"knn_inp_{p_col}",
                        )
                    )

        with st.expander(t["adv_filter_exp"], expanded=False):
            fc1, fc2, fc3 = st.columns(3)
            with fc1:
                require_complete_funnel = st.checkbox(t["req_complete"], value=True)
            with fc2:
                min_eggs_filter = st.number_input(t["min_eggs"], min_value=0, max_value=50, value=0, step=1)
            with fc3:
                use_custom_eggs = st.checkbox(t["anchor_eggs"], value=False)
                custom_eggs_val = (
                    st.number_input(t["my_eggs"], min_value=1, max_value=60, value=10, step=1)
                    if use_custom_eggs
                    else None
                )

        if knn_pool_choice == STUDY_KEY:
            candidate_pool = study_df.copy()
        elif knn_pool_choice == COMMUNITY_KEY:
            candidate_pool = comm_df.copy()
        else:
            candidate_pool = pd.concat([study_df, comm_df], ignore_index=True)

        if not candidate_pool.empty:
            if require_complete_funnel:
                candidate_pool = candidate_pool.dropna(subset=FUNNEL_STAGE_COLS)
            if min_eggs_filter > 0:
                candidate_pool = candidate_pool[candidate_pool["eggs_retrieved"] >= min_eggs_filter]

        if not selected_knn_params:
            st.warning(t["select_at_least_one_param"])
            neighbors_df = pd.DataFrame()
        else:
            neighbors_df = find_nearest_neighbors(candidate_pool, user_params=user_params, k=k_neighbors)

        if selected_knn_params and neighbors_df.empty:
            st.warning(t["no_neighbors"])
        elif not neighbors_df.empty:
            st.markdown(t["step1_neighbors"])

            base_table_cols = {
                "Include": True,
                ("排名" if lang_choice == "zh" else "Rank"): neighbors_df["rank"],
                "ID": neighbors_df["record_id"],
                ("数据集" if lang_choice == "zh" else "Dataset"): neighbors_df["dataset_key"].map(ds_names),
                ("距离" if lang_choice == "zh" else "Dist"): neighbors_df["distance"].round(2),
                p_labels["age"]: neighbors_df["age"].round(1),
                p_labels["amh"]: neighbors_df["amh"].round(2),
            }
            default_shown_keys = {
                "age",
                "amh",
                "eggs_retrieved",
                "eggs_mature",
                "fertilized",
                "blastocysts",
                "euploid",
                "mature_rate",
                "fertilization_rate",
                "blastocyst_rate",
                "euploid_rate",
                "gonadotropin_total_iu",
            }
            for extra_p in selected_knn_params:
                if extra_p not in default_shown_keys and extra_p in neighbors_df.columns:
                    base_table_cols[p_labels.get(extra_p, extra_p)] = neighbors_df[extra_p].round(2)

            base_table_cols.update(
                {
                    p_labels["eggs_retrieved"]: neighbors_df["eggs_retrieved"],
                    p_labels["eggs_mature"]: neighbors_df["eggs_mature"],
                    p_labels["fertilized"]: neighbors_df["fertilized"],
                    p_labels["blastocysts"]: neighbors_df["blastocysts"],
                    p_labels["euploid"]: neighbors_df["euploid"],
                    f"{p_labels['mature_rate']} %": (neighbors_df["mature_rate"] * 100.0).round(1),
                    f"{p_labels['fertilization_rate']} %": (neighbors_df["fertilization_rate"] * 100.0).round(1),
                    f"{p_labels['blastocyst_rate']} %": (neighbors_df["blastocyst_rate"] * 100.0).round(1),
                    f"{p_labels['euploid_rate']} %": (neighbors_df["euploid_rate"] * 100.0).round(1),
                    p_labels["gonadotropin_total_iu"]: neighbors_df["gonadotropin_total_iu"],
                    ("原因" if lang_choice == "zh" else "Reason"): neighbors_df["ivf_reason"].map(
                        lambda v: translate_ivf_reason(v, lang_choice)
                    ),
                    ("备注" if lang_choice == "zh" else "Note"): neighbors_df["note"],
                }
            )
            editor_df = pd.DataFrame(base_table_cols)

            edited_neighbors = st.data_editor(
                editor_df,
                hide_index=True,
                use_container_width=True,
                disabled=[c for c in editor_df.columns if c != "Include"],
                key="knn_neighbor_editor",
            )

            kept_ids = set(edited_neighbors.loc[edited_neighbors["Include"], "ID"])
            active_neighbors = neighbors_df[neighbors_df["record_id"].isin(kept_ids)].copy()

            st.caption(
                t["active_count"].format(
                    active=len(active_neighbors),
                    total=len(neighbors_df),
                    removed=len(neighbors_df) - len(active_neighbors),
                )
            )

            if active_neighbors.empty:
                st.error(t["all_unchecked"])
            else:
                funnel_summary = summarize_funnel(
                    active_neighbors, lang=lang_choice, starting_eggs_override=custom_eggs_val
                )

                st.markdown(t["step2_funnel"])
                kpi_cols = st.columns(len(FUNNEL_STAGE_COLS))
                for idx, stage_col in enumerate(FUNNEL_STAGE_COLS):
                    row_s = funnel_summary[funnel_summary["stage_col"] == stage_col].iloc[0]
                    mean_c = row_s["Mean"]
                    med_c = row_s["Median"]
                    cum_pct = row_s["Cumul. %"]
                    kpi_cols[idx].metric(
                        label=p_labels[stage_col],
                        value=f"{mean_c:.1f} (Med {med_c:.1f})" if pd.notna(mean_c) else "N/A",
                        delta=f"{cum_pct:.1f}{t['of_retrieved']}" if pd.notna(cum_pct) else None,
                        delta_color="off",
                    )

                viz_c1, viz_c2 = st.columns([1.2, 1.2])
                with viz_c1:
                    fig_funnel = go.Figure()
                    stage_labels = [p_labels[sc] for sc in FUNNEL_STAGE_COLS]

                    for _, nrow in active_neighbors.iterrows():
                        d_key = nrow["dataset_key"]
                        c_hex = DATASET_COLORS.get(d_key, "#555555")
                        y_vals = [nrow[c] for c in FUNNEL_STAGE_COLS]
                        age_txt = f"{nrow['age']:.1f}" if pd.notna(nrow["age"]) else "—"
                        amh_txt = f"{nrow['amh']:.2f}" if pd.notna(nrow["amh"]) else "—"
                        fig_funnel.add_trace(
                            go.Scatter(
                                x=stage_labels,
                                y=y_vals,
                                mode="lines+markers",
                                name=f"{nrow['record_id']}",
                                line=dict(color=c_hex, width=1.1),
                                marker=dict(color=c_hex, size=5),
                                opacity=0.35,
                                showlegend=False,
                                hovertemplate=(
                                    f"<b>{nrow['record_id']}</b> ({ds_names[d_key]})<br>"
                                    f"{p_labels['age']}: {age_txt}, AMH: {amh_txt}<br>"
                                    "%{x}: %{y}<extra></extra>"
                                ),
                            )
                        )

                    fig_funnel.add_trace(
                        go.Scatter(
                            x=stage_labels,
                            y=funnel_summary["Mean"],
                            mode="lines+markers+text",
                            name=t["mean_line"],
                            line=dict(color="#1f77b4", width=3.5),
                            marker=dict(size=9, color="#1f77b4"),
                            text=[f"{v:.1f}" if pd.notna(v) else "" for v in funnel_summary["Mean"]],
                            textposition="top center",
                        )
                    )
                    fig_funnel.add_trace(
                        go.Scatter(
                            x=stage_labels,
                            y=funnel_summary["Median"],
                            mode="lines+markers",
                            name=t["median_line"],
                            line=dict(color="#2ca02c", width=3, dash="dot"),
                            marker=dict(size=8, color="#2ca02c"),
                        )
                    )

                    user_stage_x = [p_labels[sc] for sc in FUNNEL_STAGE_COLS if sc in user_params]
                    user_stage_y = [user_params[sc] for sc in FUNNEL_STAGE_COLS if sc in user_params]
                    if user_stage_x:
                        fig_funnel.add_trace(
                            go.Scatter(
                                x=user_stage_x,
                                y=user_stage_y,
                                mode="markers",
                                name=t["your_input"],
                                marker=dict(
                                    color="#E60026",
                                    size=14,
                                    symbol="star",
                                    line=dict(width=1.2, color="#FFFFFF"),
                                ),
                                hovertemplate=f"<b>{t['your_input']}</b><br>%{{x}}: %{{y}}<extra></extra>",
                            )
                        )

                    fig_funnel.update_layout(
                        title=t["funnel_chart_title"],
                        template="plotly_white",
                        height=420,
                        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
                    )
                    st.plotly_chart(fig_funnel, use_container_width=True)

                with viz_c2:
                    if len(selected_knn_params) >= 2:
                        if len(selected_knn_params) > 2:
                            mc1, mc2 = st.columns(2)
                            with mc1:
                                map_x = st.selectbox(
                                    t["map_x_axis"],
                                    options=selected_knn_params,
                                    index=0,
                                    format_func=lambda k: p_labels.get(k, k),
                                    key="knn_map_x",
                                )
                            with mc2:
                                default_my = 1 if selected_knn_params[1] != map_x else 0
                                map_y = st.selectbox(
                                    t["map_y_axis"],
                                    options=selected_knn_params,
                                    index=default_my,
                                    format_func=lambda k: p_labels.get(k, k),
                                    key="knn_map_y",
                                )
                        else:
                            map_x, map_y = selected_knn_params[0], selected_knn_params[1]
                    else:
                        map_x = selected_knn_params[0]
                        map_y = "euploid" if map_x != "euploid" else "blastocysts"

                    x_label = p_labels.get(map_x, map_x)
                    y_label = p_labels.get(map_y, map_y)

                    fig_map = go.Figure()
                    for d_key in [STUDY_KEY, COMMUNITY_KEY]:
                        sub_pool = candidate_pool[candidate_pool["dataset_key"] == d_key].dropna(subset=[map_x, map_y])
                        if sub_pool.empty:
                            continue
                        mstyle = DATASET_MARKER_STYLE[d_key]
                        fig_map.add_trace(
                            go.Scatter(
                                x=sub_pool[map_x],
                                y=sub_pool[map_y],
                                mode="markers",
                                name=ds_names[d_key],
                                marker=dict(
                                    color=DATASET_COLORS[d_key],
                                    size=mstyle["size"],
                                    opacity=mstyle["opacity"],
                                ),
                                hoverinfo="skip",
                            )
                        )

                    for d_key in [STUDY_KEY, COMMUNITY_KEY]:
                        sub_nb = active_neighbors[active_neighbors["dataset_key"] == d_key].dropna(subset=[map_x, map_y])
                        if sub_nb.empty:
                            continue
                        fig_map.add_trace(
                            go.Scatter(
                                x=sub_nb[map_x],
                                y=sub_nb[map_y],
                                mode="markers",
                                name=f"{t['selected_nb']}: {ds_names[d_key]}",
                                marker=dict(
                                    color=DATASET_COLORS[d_key],
                                    size=11,
                                    opacity=0.95,
                                    line=dict(width=2, color="#00E5FF"),
                                ),
                                text=sub_nb["record_id"],
                                hovertemplate=(
                                    "<b>%{text}</b> (" + ds_names[d_key] + ")<br>"
                                    f"{x_label}: %{{x:.3g}} | {y_label}: %{{y:.3g}}<extra></extra>"
                                ),
                            )
                        )

                    if map_x in user_params and map_y in user_params:
                        ux_val = user_params[map_x]
                        uy_val = user_params[map_y]
                        fig_map.add_trace(
                            go.Scatter(
                                x=[ux_val],
                                y=[uy_val],
                                mode="markers",
                                name=t["your_input"],
                                marker=dict(
                                    color="#E60026",
                                    size=15,
                                    symbol="star",
                                    line=dict(width=1.2, color="#FFFFFF"),
                                ),
                                hovertemplate=(
                                    f"<b>{t['your_input']}</b><br>"
                                    f"{x_label}: {ux_val:.3g}<br>{y_label}: {uy_val:.3g}<extra></extra>"
                                ),
                            )
                        )
                    elif map_x in user_params:
                        fig_map.add_vline(
                            x=user_params[map_x],
                            line_dash="dash",
                            line_color="#E60026",
                            line_width=2,
                            annotation_text=f"{t['your_input']}: {user_params[map_x]:.3g}",
                            annotation_position="top right",
                        )

                    fig_map.update_layout(
                        title=t["map_chart_title"].format(x_label=x_label, y_label=y_label),
                        xaxis_title=x_label,
                        yaxis_title=y_label,
                        template="plotly_white",
                        height=420,
                        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
                    )
                    st.plotly_chart(fig_map, use_container_width=True)

                st.markdown(t["funnel_table_title"])
                funnel_col_map = (
                    {
                        "Stage": "阶段",
                        "N": "样本量 N",
                        "Mean": "均值",
                        "Median": "中位数",
                        "Q25": "25%分位",
                        "Q75": "75%分位",
                        "Min": "最小值",
                        "Max": "最大值",
                        "Step % (Mean)": "单步转化率(均值)",
                        "Step % (Median)": "单步转化率(中位)",
                        "Cumul. %": "累计占取卵%",
                        "Scaled Est.": "折算预测值",
                    }
                    if lang_choice == "zh"
                    else {c: c for c in funnel_summary.columns}
                )
                funnel_disp = funnel_summary.drop(columns=["stage_col"]).rename(columns=funnel_col_map)
                st.dataframe(
                    funnel_disp.style.format(
                        {
                            funnel_col_map["Mean"]: "{:.1f}",
                            funnel_col_map["Median"]: "{:.1f}",
                            funnel_col_map["Q25"]: "{:.1f}",
                            funnel_col_map["Q75"]: "{:.1f}",
                            funnel_col_map["Min"]: "{:.0f}",
                            funnel_col_map["Max"]: "{:.0f}",
                            funnel_col_map["Step % (Mean)"]: "{:.1f}%",
                            funnel_col_map["Step % (Median)"]: "{:.1f}%",
                            funnel_col_map["Cumul. %"]: "{:.1f}%",
                            **(
                                {funnel_col_map["Scaled Est."]: "{:.1f}"}
                                if "Scaled Est." in funnel_summary.columns
                                else {}
                            ),
                        },
                        na_rep="—",
                    ),
                    use_container_width=True,
                    hide_index=True,
                )

    # =========================================================================
    # TAB 2: Correlation Analysis & Interactive Scatterplot
    # =========================================================================
    with tab_corr:
        st.subheader(t["corr_header"])

        available_params = [
            p
            for p in p_labels.keys()
            if (not active_df.empty and p in active_df.columns and active_df[p].notna().sum() >= 5)
        ]
        if len(available_params) < 2:
            available_params = ["age", "amh", "eggs_retrieved", "blastocysts", "euploid", "euploid_rate"]

        col_c1, col_c2, col_c3 = st.columns([1.1, 1.1, 1.5])
        with col_c1:
            default_x = available_params.index("amh") if "amh" in available_params else 0
            param_x = st.selectbox(
                t["x_param"],
                options=available_params,
                index=default_x,
                format_func=lambda k: p_labels.get(k, k),
            )
        with col_c2:
            default_y = (
                available_params.index("euploid_rate")
                if "euploid_rate" in available_params
                else min(1, len(available_params) - 1)
            )
            param_y = st.selectbox(
                t["y_param"],
                options=available_params,
                index=default_y,
                format_func=lambda k: p_labels.get(k, k),
            )
        with col_c3:
            stats_ds_choice = st.radio(
                t["stats_ds"],
                options=[STUDY_KEY, COMMUNITY_KEY, "both"],
                format_func=lambda k: (
                    ds_names[STUDY_KEY]
                    if k == STUDY_KEY
                    else (ds_names[COMMUNITY_KEY] if k == COMMUNITY_KEY else t["compare_both"])
                ),
                index=0,
                horizontal=True,
            )

        def _render_corr_metrics(ds_key_label: str, corr_table: pd.DataFrame, px_col: str, py_col: str) -> None:
            st.markdown(f"**📊 {ds_key_label}**")
            if px_col == py_col:
                st.info(t["same_param_info"])
                return
            row_stats = lookup_pair_stats(corr_table, px_col, py_col)
            if row_stats is None:
                st.warning(t["insufficient_data"])
                return

            m1, m2, m3, m4, m5 = st.columns(5)
            r_val = row_stats["pearson_r"]
            p_raw = row_stats["pearson_p"]
            q_bh = row_stats["pearson_q_bh"]
            rho_val = row_stats["spearman_rho"]
            q_sp_bh = row_stats["spearman_q_bh"]
            n_val = int(row_stats["n"])

            pearson_sig = sig_stars(q_bh)
            spearman_sig = sig_stars(q_sp_bh)

            m1.metric(t["metric_n"], f"{n_val:,}")
            m2.metric(t["metric_r"], f"{r_val:+.2f} ({pearson_sig})")
            m3.metric(t["metric_p"], f"{p_raw:.2f} ({sig_stars(p_raw)})")
            m4.metric(t["metric_q"], f"{q_bh:.2f} ({pearson_sig})")
            m5.metric(t["metric_rho"], f"{rho_val:+.2f} ({spearman_sig})")

        if stats_ds_choice == STUDY_KEY:
            _render_corr_metrics(ds_names[STUDY_KEY], study_corr_df, param_x, param_y)
        elif stats_ds_choice == COMMUNITY_KEY:
            _render_corr_metrics(ds_names[COMMUNITY_KEY], comm_corr_df, param_x, param_y)
        else:
            c_left, c_right = st.columns(2)
            with c_left:
                _render_corr_metrics(ds_names[STUDY_KEY], study_corr_df, param_x, param_y)
            with c_right:
                _render_corr_metrics(ds_names[COMMUNITY_KEY], comm_corr_df, param_x, param_y)

        # Standard * Significance Reference Table
        st.caption(t["sig_ref_title"])

        show_trendline = st.checkbox(t["show_trend"], value=True)

        if not active_df.empty and param_x in active_df.columns and param_y in active_df.columns:
            plot_df = active_df.dropna(subset=[param_x, param_y]).copy()
            if not plot_df.empty:
                fig_scatter = go.Figure()
                for ds_key in [STUDY_KEY, COMMUNITY_KEY]:
                    sub_ds = plot_df[plot_df["dataset_key"] == ds_key]
                    if sub_ds.empty:
                        continue
                    color_hex = DATASET_COLORS[ds_key]
                    mstyle = DATASET_MARKER_STYLE[ds_key]
                    ds_label = ds_names[ds_key]

                    fig_scatter.add_trace(
                        go.Scatter(
                            x=sub_ds[param_x],
                            y=sub_ds[param_y],
                            mode="markers",
                            name=f"{ds_label} (n={len(sub_ds)})",
                            marker=dict(
                                color=color_hex,
                                size=mstyle["size"],
                                opacity=mstyle["opacity"],
                                line=dict(width=mstyle["line_width"], color="#FFFFFF"),
                            ),
                            customdata=np.stack(
                                [
                                    sub_ds["record_id"].fillna(""),
                                    sub_ds["age"].fillna(np.nan),
                                    sub_ds["amh"].fillna(np.nan),
                                    sub_ds["eggs_retrieved"].fillna(np.nan),
                                    sub_ds["blastocysts"].fillna(np.nan),
                                    sub_ds["euploid"].fillna(np.nan),
                                ],
                                axis=-1,
                            ),
                            hovertemplate=(
                                "<b>%{customdata[0]}</b> (" + ds_label + ")<br>"
                                + f"{p_labels.get(param_x, param_x)}: %{{x:.3g}}<br>"
                                + f"{p_labels.get(param_y, param_y)}: %{{y:.3g}}<br>"
                                + f"{p_labels['age']}: %{{customdata[1]:.1f}} | AMH: %{{customdata[2]:.2f}}<br>"
                                + f"{p_labels['eggs_retrieved']}: %{{customdata[3]:.0f}} → "
                                + f"{p_labels['blastocysts']}: %{{customdata[4]:.0f}} → "
                                + f"{p_labels['euploid']}: %{{customdata[5]:.0f}}<extra></extra>"
                            ),
                        )
                    )

                    if show_trendline and len(sub_ds) >= 3 and sub_ds[param_x].nunique() > 1:
                        slope, intercept, _, _, _ = stats.linregress(
                            sub_ds[param_x].to_numpy(dtype=float),
                            sub_ds[param_y].to_numpy(dtype=float),
                        )
                        x_line = np.linspace(sub_ds[param_x].min(), sub_ds[param_x].max(), 100)
                        y_line = slope * x_line + intercept
                        fig_scatter.add_trace(
                            go.Scatter(
                                x=x_line,
                                y=y_line,
                                mode="lines",
                                name=f"{ds_label} Trend",
                                line=dict(color=color_hex, width=2, dash="dash"),
                                hoverinfo="skip",
                            )
                        )

                fig_scatter.update_layout(
                    title=f"{p_labels.get(param_y, param_y)} vs. {p_labels.get(param_x, param_x)}",
                    xaxis_title=p_labels.get(param_x, param_x),
                    yaxis_title=p_labels.get(param_y, param_y),
                    template="plotly_white",
                    height=500,
                    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
                )
                st.plotly_chart(fig_scatter, use_container_width=True)

        # Full All-Pairs Correlation Table
        st.markdown("---")
        st.subheader(t["all_pairs_header"])
        target_corr_key = COMMUNITY_KEY if stats_ds_choice == COMMUNITY_KEY else STUDY_KEY
        selected_corr_table = comm_corr_df if target_corr_key == COMMUNITY_KEY else study_corr_df

        if selected_corr_table is not None and not selected_corr_table.empty:
            view_corr = selected_corr_table.copy()
            view_corr["label_x"] = view_corr["param_x"].map(lambda c: p_labels.get(c, c))
            view_corr["label_y"] = view_corr["param_y"].map(lambda c: p_labels.get(c, c))

            fc1, fc2 = st.columns([1.5, 1.5])
            with fc1:
                all_p_opts = [t["all_params"]] + sorted(
                    list(set(view_corr["label_x"]).union(set(view_corr["label_y"])))
                )
                filter_param = st.selectbox(t["filter_param"], options=all_p_opts, index=0)
            with fc2:
                only_sig = st.checkbox(t["only_sig"], value=False)

            if filter_param != t["all_params"]:
                view_corr = view_corr[
                    (view_corr["label_x"] == filter_param) | (view_corr["label_y"] == filter_param)
                ]
            if only_sig:
                view_corr = view_corr[view_corr["pearson_q_bh"] < 0.05]

            view_corr["sig"] = view_corr["pearson_q_bh"].map(sig_stars)

            display_corr = view_corr[
                [
                    "label_x",
                    "label_y",
                    "n",
                    "pearson_r",
                    "pearson_p",
                    "pearson_q_bh",
                    "sig",
                    "spearman_rho",
                    "spearman_p",
                    "spearman_q_bh",
                ]
            ].rename(
                columns={
                    "label_x": t["col_p1"],
                    "label_y": t["col_p2"],
                    "n": "N",
                    "pearson_r": "Pearson r",
                    "pearson_p": "Pearson p",
                    "pearson_q_bh": "Pearson BH q",
                    "sig": t["col_sig"],
                    "spearman_rho": "Spearman ρ",
                    "spearman_p": "Spearman p",
                    "spearman_q_bh": "Spearman BH q",
                }
            )
            st.dataframe(
                display_corr.style.format(
                    {
                        "Pearson r": "{:+.2f}",
                        "Pearson p": "{:.2f}",
                        "Pearson BH q": "{:.2f}",
                        "Spearman ρ": "{:+.2f}",
                        "Spearman p": "{:.2f}",
                        "Spearman BH q": "{:.2f}",
                    }
                ),
                use_container_width=True,
                height=340,
                hide_index=True,
            )

            # Standard * Significance Reference Table
            sig_ref_df = pd.DataFrame(
                {
                    ("符号" if lang_choice == "zh" else "Symbol"): ["****", "***", "**", "*", "ns"],
                    ("判定标准 (p / BH q)" if lang_choice == "zh" else "Criterion (p / BH q)"): [
                        "< 0.0001",
                        "< 0.001",
                        "< 0.01",
                        "< 0.05",
                        "≥ 0.05",
                    ],
                    ("含义" if lang_choice == "zh" else "Meaning"): (
                        ["极显著", "高度显著", "非常显著", "显著", "不显著"]
                        if lang_choice == "zh"
                        else [
                            "Extremely significant",
                            "Highly significant",
                            "Very significant",
                            "Significant",
                            "Not significant",
                        ]
                    ),
                }
            )
            st.table(sig_ref_df.set_index("符号" if lang_choice == "zh" else "Symbol"))

    # =========================================================================
    # TAB 3: Summary Statistics by Dataset & Raw Harmonized Data
    # =========================================================================
    with tab_data:
        all_df = pd.concat([study_df, comm_df], ignore_index=True)
        if all_df.empty:
            st.info(t["insufficient_data"])
        else:
            st.markdown(t["summary_header"])

            sc1, sc2 = st.columns(2)
            with sc1:
                sum_ds_choice = st.radio(
                    t["summary_ds_select"],
                    options=[STUDY_KEY, COMMUNITY_KEY, "both"],
                    format_func=lambda k: (
                        ds_names[STUDY_KEY]
                        if k == STUDY_KEY
                        else (ds_names[COMMUNITY_KEY] if k == COMMUNITY_KEY else t["both_name"])
                    ),
                    index=2,
                    horizontal=True,
                )
            with sc2:
                stat_label_map = (
                    {"mean": "均值 (Mean)", "median": "中位数 (Median)", "std": "标准差 (SD)", "all": "全部 (All)"}
                    if lang_choice == "zh"
                    else {"mean": "Mean", "median": "Median", "std": "SD", "all": "All"}
                )
                sum_stat_choice = st.radio(
                    t["summary_stat_select"],
                    options=["mean", "median", "std", "all"],
                    format_func=lambda k: stat_label_map[k],
                    index=0,
                    horizontal=True,
                )

            selected_ds_keys = (
                [STUDY_KEY, COMMUNITY_KEY]
                if sum_ds_choice == "both"
                else [sum_ds_choice]
            )
            selected_stats = (
                ["mean", "median", "std"]
                if sum_stat_choice == "all"
                else [sum_stat_choice]
            )

            wide_desc = build_wide_summary_table(
                all_df,
                lang=lang_choice,
                ds_names=ds_names,
                selected_ds_keys=selected_ds_keys,
                selected_stats=selected_stats,
            )
            st.dataframe(wide_desc.style.format("{:.1f}", na_rep="—"), use_container_width=True)

            st.markdown(t["raw_header"])
            display_raw = all_df[all_df["dataset_key"].isin(selected_ds_keys)].copy()
            display_raw["dataset"] = display_raw["dataset_key"].map(ds_names)
            display_raw["ivf_reason"] = display_raw["ivf_reason"].map(
                lambda v: translate_ivf_reason(v, lang_choice)
            )
            display_raw["protocol"] = display_raw["protocol"].map(
                lambda v: translate_protocol(v, lang_choice)
            )
            display_raw = display_raw.drop(columns=["dataset_key"])
            raw_col_map = {
                **p_labels,
                "record_id": "ID",
                "dataset": "数据集" if lang_choice == "zh" else "Dataset",
                "ivf_reason": "原因" if lang_choice == "zh" else "Reason",
                "protocol": "方案" if lang_choice == "zh" else "Protocol",
                "trigger": "夜针" if lang_choice == "zh" else "Trigger",
                "note": "备注" if lang_choice == "zh" else "Note",
            }
            display_raw = display_raw.rename(columns=raw_col_map)
            num_cols_raw = display_raw.select_dtypes(include=[np.number]).columns
            st.dataframe(
                display_raw.style.format(
                    {col: "{:.1f}" for col in num_cols_raw},
                    na_rep="—",
                ),
                use_container_width=True,
                height=420,
            )


if __name__ == "__main__":
    main()
