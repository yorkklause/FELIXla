#!/usr/bin/env python3
# designed by Kai, implemented by codex

from pathlib import Path


APP = Path("/app")


def read(path):
    return path.read_text()


def write(path, text):
    path.write_text(text)


def must_replace(text, old, new, label, count=1):
    found = text.count(old)
    if found < count:
        raise RuntimeError(f"{label}: expected at least {count} occurrence(s), found {found}")
    return text.replace(old, new, count)


def replace_once_if_missing(text, old, new, needle, label):
    if needle in text:
        return text
    return must_replace(text, old, new, label, 1)


def replace_between(text, start, end, new, label):
    start_pos = text.find(start)
    if start_pos < 0:
        raise RuntimeError(f"{label}: start marker not found")
    end_pos = text.find(end, start_pos)
    if end_pos < 0:
        raise RuntimeError(f"{label}: end marker not found")
    return text[:start_pos] + new + text[end_pos:]


def patch_main_hpp():
    path = APP / "src" / "Main.hpp"
    text = read(path)
    text = replace_once_if_missing(
        text,
        "bool check_Vcf_end();\n",
        "bool check_Vcf_end();\n\nbool check_TRACTORHYBRID_end();\n",
        "check_TRACTORHYBRID_end",
        "Main.hpp check_TRACTORHYBRID_end prototype",
    )
    write(path, text)


def patch_makevars():
    path = APP / "src" / "Makevars"
    text = read(path)
    if "TractorHybrid.o" not in text:
        text = must_replace(
            text,
            "Duplicate_vec.o ",
            "Duplicate_vec.o TractorHybrid.o ",
            "Makevars TractorHybrid object",
            1,
        )
    write(path, text)


def patch_main_cpp():
    path = APP / "src" / "Main.cpp"
    text = read(path)

    text = replace_once_if_missing(
        text,
        '#include "VCF.hpp"\n',
        '#include "VCF.hpp"\n#include "TractorHybrid.hpp"\n',
        'include "TractorHybrid.hpp"',
        "Main.cpp include",
    )

    text = replace_once_if_missing(
        text,
        "static VCF::VcfClass* ptr_gVCFobj = NULL;\n",
        """static VCF::VcfClass* ptr_gVCFobj = NULL;
static TRACTOR_HYBRID::TractorHybridClass* ptr_gTRACTORHYBRIDobj = NULL;
static double g_TRACTORHYBRID_get_marker_seconds = 0.0;
static double g_TRACTORHYBRID_marker_pvalue_seconds = 0.0;
static double g_TRACTORHYBRID_output_seconds = 0.0;
static double g_TRACTORHYBRID_reset_seconds = 0.0;
static double g_TRACTORHYBRID_condition_seconds = 0.0;
static double g_TRACTORHYBRID_impute_qc_seconds = 0.0;
static double g_TRACTORHYBRID_joint_seconds = 0.0;
static double g_TRACTORHYBRID_variance_seconds = 0.0;
static double g_TRACTORHYBRID_condition_cache_hits = 0.0;
static double g_TRACTORHYBRID_condition_cache_misses = 0.0;
static bool g_TRACTORHYBRID_condition_cache_valid = false;
static uint32_t g_TRACTORHYBRID_condition_cache_block_id = UINT32_MAX;
static int g_TRACTORHYBRID_condition_cache_n_anc = -1;
static double g_TRACTORHYBRID_condition_cache_cutoff = 0.0;
static std::string g_TRACTORHYBRID_condition_cache_trait;
static arma::vec g_TRACTORHYBRID_condition_cache_nanc_case_vec;
static arma::vec g_TRACTORHYBRID_condition_cache_nanc_ctrl_vec;
static arma::vec g_TRACTORHYBRID_condition_cache_nanc_vec;
static arma::uvec g_TRACTORHYBRID_condition_cache_not_nan_anc_indices_vec;
static bool g_TRACTORHYBRID_condition_cache_isconditiononHaplo = false;
static arma::mat g_TRACTORHYBRID_condition_cache_P2Mat_cond;
static arma::mat g_TRACTORHYBRID_condition_cache_VarInvMat_cond;
static arma::mat g_TRACTORHYBRID_condition_cache_VarMat_cond;
static arma::vec g_TRACTORHYBRID_condition_cache_Tstat_cond;
static arma::vec g_TRACTORHYBRID_condition_cache_G2_Weight_cond;
static arma::vec g_TRACTORHYBRID_condition_cache_MAF_cond;
static double g_TRACTORHYBRID_condition_cache_qsum_cond = 0.0;
static arma::vec g_TRACTORHYBRID_condition_cache_gsum_cond;
static std::vector<std::string> g_TRACTORHYBRID_condition_cache_p_cond;

static inline void tractor_hybrid_clear_condition_cache(){
  g_TRACTORHYBRID_condition_cache_valid = false;
  g_TRACTORHYBRID_condition_cache_block_id = UINT32_MAX;
  g_TRACTORHYBRID_condition_cache_n_anc = -1;
  g_TRACTORHYBRID_condition_cache_trait.clear();
}

static inline double tractor_hybrid_now_seconds(){
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

#define TRACTORHYBRID_TIMED_GET_MARKER_PVAL(...) do { \
  double __tractor_pvalue_start = tractor_hybrid_now_seconds(); \
  (*ptr_gSAIGEobj).getMarkerPval(__VA_ARGS__); \
  g_TRACTORHYBRID_marker_pvalue_seconds += tractor_hybrid_now_seconds() - __tractor_pvalue_start; \
} while(0)

#define TRACTORHYBRID_TIMED_ASSIGN_CONDITION(...) do { \
  double __tractor_condition_start = tractor_hybrid_now_seconds(); \
  assign_conditionHaplotypes(__VA_ARGS__); \
  g_TRACTORHYBRID_condition_seconds += tractor_hybrid_now_seconds() - __tractor_condition_start; \
} while(0)

#define TRACTORHYBRID_TIMED_IMPUTE_FAKEFLIP(...) ([&](){ \
  double __tractor_impute_start = tractor_hybrid_now_seconds(); \
  auto __tractor_impute_ret = (imputeGenoAndFlip_fakeflip)(__VA_ARGS__); \
  g_TRACTORHYBRID_impute_qc_seconds += tractor_hybrid_now_seconds() - __tractor_impute_start; \
  return __tractor_impute_ret; \
}())

#define TRACTORHYBRID_TIMED_IMPUTE(...) ([&](){ \
  double __tractor_impute_start = tractor_hybrid_now_seconds(); \
  auto __tractor_impute_ret = (imputeGenoAndFlip)(__VA_ARGS__); \
  g_TRACTORHYBRID_impute_qc_seconds += tractor_hybrid_now_seconds() - __tractor_impute_start; \
  return __tractor_impute_ret; \
}())

#define TRACTORHYBRID_TIMED_GET_JOINT(...) ([&](){ \
  double __tractor_joint_start = tractor_hybrid_now_seconds(); \
  auto __tractor_joint_ret = (get_jointScore_pvalue)(__VA_ARGS__); \
  g_TRACTORHYBRID_joint_seconds += tractor_hybrid_now_seconds() - __tractor_joint_start; \
  return __tractor_joint_ret; \
}())

#define TRACTORHYBRID_TIMED_CCT(...) ([&](){ \
  double __tractor_joint_start = tractor_hybrid_now_seconds(); \
  auto __tractor_joint_ret = (CCT_cpp)(__VA_ARGS__); \
  g_TRACTORHYBRID_joint_seconds += tractor_hybrid_now_seconds() - __tractor_joint_start; \
  return __tractor_joint_ret; \
}())

#define TRACTORHYBRID_TIMED_SPA_ER(...) do { \
  double __tractor_joint_start = tractor_hybrid_now_seconds(); \
  (SPA_ER_kernel_related_Phiadj_fast_new_cpp)(__VA_ARGS__); \
  g_TRACTORHYBRID_joint_seconds += tractor_hybrid_now_seconds() - __tractor_joint_start; \
} while(0)

#define TRACTORHYBRID_TIMED_SET_SPARSE_FLAG(...) do { \
  double __tractor_variance_start = tractor_hybrid_now_seconds(); \
  (*ptr_gSAIGEobj).set_flagSparseGRM_cur(__VA_ARGS__); \
  g_TRACTORHYBRID_variance_seconds += tractor_hybrid_now_seconds() - __tractor_variance_start; \
} while(0)

#define TRACTORHYBRID_TIMED_ASSIGN_SINGLE_VARIANCE(...) do { \
  double __tractor_variance_start = tractor_hybrid_now_seconds(); \
  (*ptr_gSAIGEobj).assignSingleVarianceRatio(__VA_ARGS__); \
  g_TRACTORHYBRID_variance_seconds += tractor_hybrid_now_seconds() - __tractor_variance_start; \
} while(0)

#define TRACTORHYBRID_TIMED_ASSIGN_VARIANCE(...) ([&](){ \
  double __tractor_variance_start = tractor_hybrid_now_seconds(); \
  auto __tractor_variance_ret = (*ptr_gSAIGEobj).assignVarianceRatio(__VA_ARGS__); \
  g_TRACTORHYBRID_variance_seconds += tractor_hybrid_now_seconds() - __tractor_variance_start; \
  return __tractor_variance_ret; \
}())
""",
        "ptr_gTRACTORHYBRIDobj",
        "Main.cpp global pointer",
    )

    text = replace_once_if_missing(
        text,
        "  if(g_is_rewrite_XnonPAR_forMales){\n",
        """  if(t_genoType == "tractor_hybrid"){
    double __tractor_get_marker_start = tractor_hybrid_now_seconds();
    ptr_gTRACTORHYBRIDobj->getOneMarker(t_gIndex_prev, t_gIndex, t_ref, t_alt, t_marker, t_pd, t_chr,
                                      t_altFreq, t_altCounts, t_missingRate, t_imputeInfo,
                                      t_isOutputIndexForMissing, t_indexForMissing,
                                      t_isOnlyOutputNonZero, t_indexForNonZero,
                                      isBoolRead, t_GVec, t_isImputation);
    g_TRACTORHYBRID_get_marker_seconds += tractor_hybrid_now_seconds() - __tractor_get_marker_start;
  }

  if(g_is_rewrite_XnonPAR_forMales){
""",
        't_genoType == "tractor_hybrid"',
        "Main.cpp Unified_getOneMarker branch",
    )

    text = replace_once_if_missing(
        text,
        """    if(t_genoType == "vcf"){
       N0 = ptr_gVCFobj->getN0();
    }
    return(N0);
""",
        """    if(t_genoType == "vcf"){
       N0 = ptr_gVCFobj->getN0();
    }
    if(t_genoType == "tractor_hybrid"){
       N0 = ptr_gTRACTORHYBRIDobj->getN0();
    }
    return(N0);
""",
        "ptr_gTRACTORHYBRIDobj->getN0",
        "Main.cpp getN0 branch",
    )

    text = replace_once_if_missing(
        text,
        """    if(t_genoType == "vcf"){
       N = ptr_gVCFobj->getN();
    }
    return(N);
""",
        """    if(t_genoType == "vcf"){
       N = ptr_gVCFobj->getN();
    }
    if(t_genoType == "tractor_hybrid"){
       N = ptr_gTRACTORHYBRIDobj->getN();
    }
    return(N);
""",
        "ptr_gTRACTORHYBRIDobj->getN()",
        "Main.cpp getN branch",
    )

    text = replace_once_if_missing(
        text,
        """bool isEnd = ptr_gVCFobj->check_iterator_end();

}
""",
        """bool isEnd = ptr_gVCFobj->check_iterator_end();

}

// [[Rcpp::export]]
void setTRACTORHYBRIDobjInCPP(std::string t_prefix,
            std::vector<std::string> & t_SampleInModel)
{
  tractor_hybrid_clear_condition_cache();
  ptr_gTRACTORHYBRIDobj = new TRACTOR_HYBRID::TractorHybridClass(t_prefix, t_SampleInModel);
}
""",
        "setTRACTORHYBRIDobjInCPP",
        "Main.cpp hybrid constructor",
    )

    text = replace_once_if_missing(
        text,
        """  if(t_genoType == "vcf"){
    ptr_gVCFobj->m_fmtField = t_vcfField;
    //std::cout << "ptr_gVCFobj->m_fmtField " << ptr_gVCFobj->m_fmtField << std::endl;
    ptr_gVCFobj->getOneMarker(t_ref, t_alt, t_marker, t_pd, t_chr, t_altFreq, t_altCounts, t_missingRate, t_imputeInfo,
                                      t_isOutputIndexForMissing, t_indexForMissing, t_isOnlyOutputNonZero, t_indexForNonZero, isBoolRead, t_GVec, t_isImputation);
    //ptr_gVCFobj->move_forward_iterator(1);
  }

  if(g_is_rewrite_XnonPAR_forMales){
""",
        """  if(t_genoType == "vcf"){
    ptr_gVCFobj->m_fmtField = t_vcfField;
    //std::cout << "ptr_gVCFobj->m_fmtField " << ptr_gVCFobj->m_fmtField << std::endl;
    ptr_gVCFobj->getOneMarker(t_ref, t_alt, t_marker, t_pd, t_chr, t_altFreq, t_altCounts, t_missingRate, t_imputeInfo,
                                      t_isOutputIndexForMissing, t_indexForMissing, t_isOnlyOutputNonZero, t_indexForNonZero, isBoolRead, t_GVec, t_isImputation);
    //ptr_gVCFobj->move_forward_iterator(1);
  }

  if(t_genoType == "tractor_hybrid"){
    double __tractor_get_marker_start = tractor_hybrid_now_seconds();
    ptr_gTRACTORHYBRIDobj->getOneMarkerAdmixedField(t_gIndex_prev, t_gIndex, t_ref, t_alt, t_marker, t_pd, t_chr,
                                      t_altFreq, t_altCounts, t_missingRate, t_imputeInfo,
                                      t_isOutputIndexForMissing, t_indexForMissing,
                                      t_isOnlyOutputNonZero, t_indexForNonZero,
                                      isBoolRead, t_GVec, t_isImputation, t_vcfField);
    g_TRACTORHYBRID_get_marker_seconds += tractor_hybrid_now_seconds() - __tractor_get_marker_start;
  }

  if(g_is_rewrite_XnonPAR_forMales){
""",
        "getOneMarkerAdmixedField",
        "Main.cpp Unified_getOneMarker_Admixed branch",
    )

    text = replace_once_if_missing(
        text,
        """void move_forward_iterator_Vcf(int i){
	ptr_gVCFobj->move_forward_iterator(i);
}
""",
        """void move_forward_iterator_Vcf(int i){
	ptr_gVCFobj->move_forward_iterator(i);
}

// [[Rcpp::export]]
void set_iterator_inTRACTORHYBRID(std::string chrom, int beg_pd, int end_pd){
	ptr_gTRACTORHYBRIDobj->set_iterator(chrom, beg_pd, end_pd);
}

// [[Rcpp::export]]
bool check_TRACTORHYBRID_end(){
	return ptr_gTRACTORHYBRIDobj->check_iterator_end();
}

// [[Rcpp::export]]
void move_forward_iterator_TRACTORHYBRID(int i){
	ptr_gTRACTORHYBRIDobj->move_forward_iterator(i);
}

// [[Rcpp::export]]
double get_TRACTORHYBRID_io_seconds(){
	if(ptr_gTRACTORHYBRIDobj == NULL){
		return 0.0;
	}
	return ptr_gTRACTORHYBRIDobj->get_io_seconds();
}

// [[Rcpp::export]]
double get_TRACTORHYBRID_decode_seconds(){
	if(ptr_gTRACTORHYBRIDobj == NULL){
		return 0.0;
	}
	return ptr_gTRACTORHYBRIDobj->get_decode_seconds();
}

// [[Rcpp::export]]
double get_TRACTORHYBRID_get_marker_seconds(){
	return g_TRACTORHYBRID_get_marker_seconds;
}

// [[Rcpp::export]]
double get_TRACTORHYBRID_pvalue_seconds(){
	return g_TRACTORHYBRID_marker_pvalue_seconds;
}

// [[Rcpp::export]]
double get_TRACTORHYBRID_output_seconds(){
	return g_TRACTORHYBRID_output_seconds;
}

// [[Rcpp::export]]
double get_TRACTORHYBRID_reset_seconds(){
	return g_TRACTORHYBRID_reset_seconds;
}

// [[Rcpp::export]]
double get_TRACTORHYBRID_condition_seconds(){
	return g_TRACTORHYBRID_condition_seconds;
}

// [[Rcpp::export]]
double get_TRACTORHYBRID_impute_qc_seconds(){
	return g_TRACTORHYBRID_impute_qc_seconds;
}

// [[Rcpp::export]]
double get_TRACTORHYBRID_joint_seconds(){
	return g_TRACTORHYBRID_joint_seconds;
}

// [[Rcpp::export]]
double get_TRACTORHYBRID_variance_seconds(){
	return g_TRACTORHYBRID_variance_seconds;
}

// [[Rcpp::export]]
double get_TRACTORHYBRID_condition_cache_hits(){
	return g_TRACTORHYBRID_condition_cache_hits;
}

// [[Rcpp::export]]
double get_TRACTORHYBRID_condition_cache_misses(){
	return g_TRACTORHYBRID_condition_cache_misses;
}
""",
        "set_iterator_inTRACTORHYBRID",
        "Main.cpp hybrid iterator helpers",
    )

    n_get_marker_pval = text.count("ptr_gSAIGEobj->getMarkerPval(")
    if n_get_marker_pval == 0:
        raise RuntimeError("Main.cpp marker pvalue timer: no getMarkerPval calls found")
    text = text.replace("ptr_gSAIGEobj->getMarkerPval(", "TRACTORHYBRID_TIMED_GET_MARKER_PVAL(")
    text = text.replace("    assign_conditionHaplotypes(\n", "    TRACTORHYBRID_TIMED_ASSIGN_CONDITION(\n")
    text = text.replace("flip = imputeGenoAndFlip_fakeflip(", "flip = TRACTORHYBRID_TIMED_IMPUTE_FAKEFLIP(")
    text = text.replace("flip = imputeGenoAndFlip(", "flip = TRACTORHYBRID_TIMED_IMPUTE(")
    text = text.replace("= get_jointScore_pvalue(", "= TRACTORHYBRID_TIMED_GET_JOINT(")
    text = text.replace(" CCT_cpp(", " TRACTORHYBRID_TIMED_CCT(")
    text = text.replace("SPA_ER_kernel_related_Phiadj_fast_new_cpp(", "TRACTORHYBRID_TIMED_SPA_ER(")
    text = text.replace("ptr_gSAIGEobj->set_flagSparseGRM_cur(", "TRACTORHYBRID_TIMED_SET_SPARSE_FLAG(")
    text = text.replace("ptr_gSAIGEobj->assignSingleVarianceRatio(", "TRACTORHYBRID_TIMED_ASSIGN_SINGLE_VARIANCE(")
    text = text.replace("ptr_gSAIGEobj->assignVarianceRatio(", "TRACTORHYBRID_TIMED_ASSIGN_VARIANCE(")

    text = replace_once_if_missing(
        text,
        """    bool isconditiononHaplo;
    TRACTORHYBRID_TIMED_ASSIGN_CONDITION(
    		t_traitType,
                           t_genoType,     //"vcf"
                           n,
			   t_NumberofANC,	
                           gIndex,
                           gIndex_prev,
                           nanc_case_vec,
                           nanc_ctrl_vec,
                           nanc_vec,
                           not_nan_anc_indices_vec,
			   t_pvalcutoff_of_haplotype,
			   isconditiononHaplo
                           );           // sample size
""",
        """    bool isconditiononHaplo;
    uint32_t __tractor_condition_block_id = UINT32_MAX;
    bool __tractor_condition_cache_hit = false;

    if(t_genoType == "tractor_hybrid" && ptr_gTRACTORHYBRIDobj != NULL){
      __tractor_condition_block_id = ptr_gTRACTORHYBRIDobj->blockIdForMarker(gIndex);
      if(g_TRACTORHYBRID_condition_cache_valid &&
         g_TRACTORHYBRID_condition_cache_block_id == __tractor_condition_block_id &&
         g_TRACTORHYBRID_condition_cache_n_anc == t_NumberofANC &&
         g_TRACTORHYBRID_condition_cache_cutoff == t_pvalcutoff_of_haplotype &&
         g_TRACTORHYBRID_condition_cache_trait == t_traitType){
        nanc_case_vec = g_TRACTORHYBRID_condition_cache_nanc_case_vec;
        nanc_ctrl_vec = g_TRACTORHYBRID_condition_cache_nanc_ctrl_vec;
        nanc_vec = g_TRACTORHYBRID_condition_cache_nanc_vec;
        not_nan_anc_indices_vec = g_TRACTORHYBRID_condition_cache_not_nan_anc_indices_vec;
        isconditiononHaplo = g_TRACTORHYBRID_condition_cache_isconditiononHaplo;
        ptr_gSAIGEobj->m_P2Mat_cond = g_TRACTORHYBRID_condition_cache_P2Mat_cond;
        ptr_gSAIGEobj->m_VarInvMat_cond = g_TRACTORHYBRID_condition_cache_VarInvMat_cond;
        ptr_gSAIGEobj->m_VarMat_cond = g_TRACTORHYBRID_condition_cache_VarMat_cond;
        ptr_gSAIGEobj->m_Tstat_cond = g_TRACTORHYBRID_condition_cache_Tstat_cond;
        ptr_gSAIGEobj->m_G2_Weight_cond = g_TRACTORHYBRID_condition_cache_G2_Weight_cond;
        ptr_gSAIGEobj->m_MAF_cond = g_TRACTORHYBRID_condition_cache_MAF_cond;
        ptr_gSAIGEobj->m_qsum_cond = g_TRACTORHYBRID_condition_cache_qsum_cond;
        ptr_gSAIGEobj->m_gsum_cond = g_TRACTORHYBRID_condition_cache_gsum_cond;
        ptr_gSAIGEobj->m_p_cond = g_TRACTORHYBRID_condition_cache_p_cond;
        g_TRACTORHYBRID_condition_cache_hits += 1.0;
        __tractor_condition_cache_hit = true;
      }
    }

    if(!__tractor_condition_cache_hit){
      TRACTORHYBRID_TIMED_ASSIGN_CONDITION(
    		t_traitType,
                           t_genoType,     //"vcf"
                           n,
			   t_NumberofANC,	
                           gIndex,
                           gIndex_prev,
                           nanc_case_vec,
                           nanc_ctrl_vec,
                           nanc_vec,
                           not_nan_anc_indices_vec,
			   t_pvalcutoff_of_haplotype,
			   isconditiononHaplo
                           );           // sample size

      if(t_genoType == "tractor_hybrid" && __tractor_condition_block_id != UINT32_MAX){
        g_TRACTORHYBRID_condition_cache_valid = true;
        g_TRACTORHYBRID_condition_cache_block_id = __tractor_condition_block_id;
        g_TRACTORHYBRID_condition_cache_n_anc = t_NumberofANC;
        g_TRACTORHYBRID_condition_cache_cutoff = t_pvalcutoff_of_haplotype;
        g_TRACTORHYBRID_condition_cache_trait = t_traitType;
        g_TRACTORHYBRID_condition_cache_nanc_case_vec = nanc_case_vec;
        g_TRACTORHYBRID_condition_cache_nanc_ctrl_vec = nanc_ctrl_vec;
        g_TRACTORHYBRID_condition_cache_nanc_vec = nanc_vec;
        g_TRACTORHYBRID_condition_cache_not_nan_anc_indices_vec = not_nan_anc_indices_vec;
        g_TRACTORHYBRID_condition_cache_isconditiononHaplo = isconditiononHaplo;
        g_TRACTORHYBRID_condition_cache_P2Mat_cond = ptr_gSAIGEobj->m_P2Mat_cond;
        g_TRACTORHYBRID_condition_cache_VarInvMat_cond = ptr_gSAIGEobj->m_VarInvMat_cond;
        g_TRACTORHYBRID_condition_cache_VarMat_cond = ptr_gSAIGEobj->m_VarMat_cond;
        g_TRACTORHYBRID_condition_cache_Tstat_cond = ptr_gSAIGEobj->m_Tstat_cond;
        g_TRACTORHYBRID_condition_cache_G2_Weight_cond = ptr_gSAIGEobj->m_G2_Weight_cond;
        g_TRACTORHYBRID_condition_cache_MAF_cond = ptr_gSAIGEobj->m_MAF_cond;
        g_TRACTORHYBRID_condition_cache_qsum_cond = ptr_gSAIGEobj->m_qsum_cond;
        g_TRACTORHYBRID_condition_cache_gsum_cond = ptr_gSAIGEobj->m_gsum_cond;
        g_TRACTORHYBRID_condition_cache_p_cond = ptr_gSAIGEobj->m_p_cond;
        g_TRACTORHYBRID_condition_cache_misses += 1.0;
      }
    }
""",
        "__tractor_condition_cache_hit",
        "Main.cpp tractor_hybrid condition cache",
    )

    text = replace_once_if_missing(
        text,
        """    bool isSPAjoint = false;
    t_GVecHom.zeros();
    Scorevec.fill(arma::datum::nan);
    PvalvecallAnc.fill(arma::datum::nan);
    Scorevec_c.fill(arma::datum::nan);
    //Scorevec.zeros();
    //Scorevec_c.zeros();
    P1Mat.zeros();
    P2Mat.zeros();
    G1tilde_P_G2tilde_Mat.zeros();
    numancresults = 0;
""",
        """    bool isSPAjoint = false;
    double __tractor_reset_start = tractor_hybrid_now_seconds();
    t_GVecHom.zeros();
    Scorevec.fill(arma::datum::nan);
    PvalvecallAnc.fill(arma::datum::nan);
    Scorevec_c.fill(arma::datum::nan);
    //Scorevec.zeros();
    //Scorevec_c.zeros();
    P1Mat.zeros();
    P2Mat.zeros();
    G1tilde_P_G2tilde_Mat.zeros();
    g_TRACTORHYBRID_reset_seconds += tractor_hybrid_now_seconds() - __tractor_reset_start;
    numancresults = 0;
""",
        "__tractor_reset_start",
        "Main.cpp admixed reset timer",
    )

    text = replace_once_if_missing(
        text,
        """			){
  t_OutFile_singleInGroup << markerName;
""",
        """			){
  double __tractor_output_start = tractor_hybrid_now_seconds();
  t_OutFile_singleInGroup << markerName;
""",
        "__tractor_output_start",
        "Main.cpp admixed output timer start",
    )

    text = replace_once_if_missing(
        text,
        """  //output
  writeOutfile_single_admixed_new(t_isMoreOutput,
""",
        """  //output
  double __tractor_output_new_start = tractor_hybrid_now_seconds();
  writeOutfile_single_admixed_new(t_isMoreOutput,
""",
        "__tractor_output_new_start",
        "Main.cpp admixed new output timer start",
    )

    text = replace_once_if_missing(
        text,
        """  pvalHap_Vec);

}



void writeOutfile_single_admixed_new""",
        """  pvalHap_Vec);
  g_TRACTORHYBRID_output_seconds += tractor_hybrid_now_seconds() - __tractor_output_new_start;

}



void writeOutfile_single_admixed_new""",
        "g_TRACTORHYBRID_output_seconds += tractor_hybrid_now_seconds() - __tractor_output_new_start",
        "Main.cpp admixed new output timer end",
    )

    text = replace_once_if_missing(
        text,
        """  //t_OutFile_singleInGroup << "\\n";
  return(numofAnc);
}
""",
        """  //t_OutFile_singleInGroup << "\\n";
  g_TRACTORHYBRID_output_seconds += tractor_hybrid_now_seconds() - __tractor_output_start;
  return(numofAnc);
}
""",
        "g_TRACTORHYBRID_output_seconds += tractor_hybrid_now_seconds()",
        "Main.cpp admixed output timer end",
    )

    text = replace_once_if_missing(
        text,
        """  }else if(t_genoType == "vcf"){
    ptr_gVCFobj->closegenofile();
  }else if(t_genoType == "plink"){
""",
        """  }else if(t_genoType == "vcf"){
    ptr_gVCFobj->closegenofile();
  }else if(t_genoType == "tractor_hybrid"){
    ptr_gTRACTORHYBRIDobj->closegenofile();
  }else if(t_genoType == "plink"){
""",
        "ptr_gTRACTORHYBRIDobj->closegenofile",
        "Main.cpp close branch",
    )

    text = text.replace(
        "bool isEndFile = check_Vcf_end();",
        'bool isEndFile = (t_genoType == "tractor_hybrid") ? check_TRACTORHYBRID_end() : check_Vcf_end();',
    )

    text = text.replace(
        "  ptr_gVCFobj->move_forward_iterator(1);\n  //std::cout << \"g_current_vcffield_ancestry\"",
        """  if(t_genoType == "vcf"){
    ptr_gVCFobj->move_forward_iterator(1);
  }
  if(t_genoType == "tractor_hybrid"){
    ptr_gTRACTORHYBRIDobj->move_forward_iterator(1);
  }
  //std::cout << "g_current_vcffield_ancestry\"""",
    )

    text = text.replace(
        "\n    ptr_gVCFobj->move_forward_iterator(1);\n",
        """
    if(t_genoType == "vcf"){
      ptr_gVCFobj->move_forward_iterator(1);
    }
    if(t_genoType == "tractor_hybrid"){
      ptr_gTRACTORHYBRIDobj->move_forward_iterator(1);
    }
""",
    )

    write(path, text)


def patch_saige_test_cpp():
    path = APP / "src" / "SAIGE_test.cpp"
    text = read(path)

    text = replace_once_if_missing(
        text,
        """#include <boost/math/distributions/chi_squared.hpp>

namespace SAIGE {
""",
        """#include <boost/math/distributions/chi_squared.hpp>

static double g_TRACTORHYBRID_saige_score_seconds = 0.0;
static double g_TRACTORHYBRID_saige_score_fast_calls = 0.0;
static double g_TRACTORHYBRID_saige_score_slow_calls = 0.0;
static double g_TRACTORHYBRID_saige_scorefast_extract_seconds = 0.0;
static double g_TRACTORHYBRID_saige_scorefast_projection_seconds = 0.0;
static double g_TRACTORHYBRID_saige_scorefast_variance_seconds = 0.0;
static double g_TRACTORHYBRID_saige_scorefast_result_seconds = 0.0;
static double g_TRACTORHYBRID_saige_scorefast_gtilde_seconds = 0.0;
static double g_TRACTORHYBRID_saige_scorefast_gtilde_calls = 0.0;
static arma::vec* g_TRACTORHYBRID_saige_scorefast_gtilde_target = nullptr;
static bool g_TRACTORHYBRID_saige_scorefast_gtilde_ready = false;
static double g_TRACTORHYBRID_saige_alloc_seconds = 0.0;
static double g_TRACTORHYBRID_saige_getadjg_seconds = 0.0;
static double g_TRACTORHYBRID_saige_getadjg_calls = 0.0;
static double g_TRACTORHYBRID_saige_getadjg_accum_seconds = 0.0;
static double g_TRACTORHYBRID_saige_getadjg_projection_seconds = 0.0;
static double g_TRACTORHYBRID_saige_getadjg_spa_seconds = 0.0;
static double g_TRACTORHYBRID_saige_getadjg_spa_calls = 0.0;
static double g_TRACTORHYBRID_saige_getadjg_firth_seconds = 0.0;
static double g_TRACTORHYBRID_saige_getadjg_firth_calls = 0.0;
static double g_TRACTORHYBRID_saige_getadjg_condition_seconds = 0.0;
static double g_TRACTORHYBRID_saige_getadjg_condition_calls = 0.0;
static double g_TRACTORHYBRID_saige_getadjg_region_seconds = 0.0;
static double g_TRACTORHYBRID_saige_getadjg_region_calls = 0.0;
static double g_TRACTORHYBRID_saige_getadjg_other_seconds = 0.0;
static double g_TRACTORHYBRID_saige_getadjg_other_calls = 0.0;
static double g_TRACTORHYBRID_saige_spa_seconds = 0.0;
static double g_TRACTORHYBRID_saige_spa_calls = 0.0;
static double g_TRACTORHYBRID_saige_firth_seconds = 0.0;
static double g_TRACTORHYBRID_saige_firth_calls = 0.0;
static double g_TRACTORHYBRID_saige_er_seconds = 0.0;
static double g_TRACTORHYBRID_saige_condition_seconds = 0.0;
static double g_TRACTORHYBRID_saige_region_seconds = 0.0;

static inline double tractor_hybrid_saige_now_seconds(){
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

#define TRACTORHYBRID_SAIGE_TIMED_SCORE(...) do { \\
  double __tractor_saige_start = tractor_hybrid_saige_now_seconds(); \\
  scoreTest(__VA_ARGS__); \\
  g_TRACTORHYBRID_saige_score_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_start; \\
  g_TRACTORHYBRID_saige_score_slow_calls += 1.0; \\
} while(0)

#define TRACTORHYBRID_SAIGE_TIMED_SCORE_FAST(...) do { \\
  double __tractor_saige_start = tractor_hybrid_saige_now_seconds(); \\
  scoreTestFast(__VA_ARGS__); \\
  g_TRACTORHYBRID_saige_score_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_start; \\
  g_TRACTORHYBRID_saige_score_fast_calls += 1.0; \\
} while(0)

#define TRACTORHYBRID_SAIGE_TIMED_GETADJG_REASON(REASON_SECONDS, REASON_CALLS, ...) do { \\
  double __tractor_saige_start = tractor_hybrid_saige_now_seconds(); \\
  getadjGFast(__VA_ARGS__); \\
  double __tractor_saige_delta = tractor_hybrid_saige_now_seconds() - __tractor_saige_start; \\
  g_TRACTORHYBRID_saige_getadjg_seconds += __tractor_saige_delta; \\
  g_TRACTORHYBRID_saige_getadjg_calls += 1.0; \\
  REASON_SECONDS += __tractor_saige_delta; \\
  REASON_CALLS += 1.0; \\
} while(0)

#define TRACTORHYBRID_SAIGE_TIMED_GETADJG(...) \\
  TRACTORHYBRID_SAIGE_TIMED_GETADJG_REASON(g_TRACTORHYBRID_saige_getadjg_other_seconds, g_TRACTORHYBRID_saige_getadjg_other_calls, __VA_ARGS__)

#define TRACTORHYBRID_SAIGE_TIMED_SPA_FAST(...) do { \\
  double __tractor_saige_start = tractor_hybrid_saige_now_seconds(); \\
  SPA_fast(__VA_ARGS__); \\
  g_TRACTORHYBRID_saige_spa_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_start; \\
  g_TRACTORHYBRID_saige_spa_calls += 1.0; \\
} while(0)

#define TRACTORHYBRID_SAIGE_TIMED_SPA(...) do { \\
  double __tractor_saige_start = tractor_hybrid_saige_now_seconds(); \\
  SPA(__VA_ARGS__); \\
  g_TRACTORHYBRID_saige_spa_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_start; \\
  g_TRACTORHYBRID_saige_spa_calls += 1.0; \\
} while(0)

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_score_seconds(){ return g_TRACTORHYBRID_saige_score_seconds; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_score_fast_calls(){ return g_TRACTORHYBRID_saige_score_fast_calls; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_score_slow_calls(){ return g_TRACTORHYBRID_saige_score_slow_calls; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_scorefast_extract_seconds(){ return g_TRACTORHYBRID_saige_scorefast_extract_seconds; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_scorefast_projection_seconds(){ return g_TRACTORHYBRID_saige_scorefast_projection_seconds; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_scorefast_variance_seconds(){ return g_TRACTORHYBRID_saige_scorefast_variance_seconds; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_scorefast_result_seconds(){ return g_TRACTORHYBRID_saige_scorefast_result_seconds; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_scorefast_gtilde_seconds(){ return g_TRACTORHYBRID_saige_scorefast_gtilde_seconds; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_scorefast_gtilde_calls(){ return g_TRACTORHYBRID_saige_scorefast_gtilde_calls; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_alloc_seconds(){ return g_TRACTORHYBRID_saige_alloc_seconds; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_getadjg_seconds(){ return g_TRACTORHYBRID_saige_getadjg_seconds; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_getadjg_calls(){ return g_TRACTORHYBRID_saige_getadjg_calls; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_getadjg_accum_seconds(){ return g_TRACTORHYBRID_saige_getadjg_accum_seconds; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_getadjg_projection_seconds(){ return g_TRACTORHYBRID_saige_getadjg_projection_seconds; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_getadjg_spa_seconds(){ return g_TRACTORHYBRID_saige_getadjg_spa_seconds; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_getadjg_spa_calls(){ return g_TRACTORHYBRID_saige_getadjg_spa_calls; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_getadjg_firth_seconds(){ return g_TRACTORHYBRID_saige_getadjg_firth_seconds; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_getadjg_firth_calls(){ return g_TRACTORHYBRID_saige_getadjg_firth_calls; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_getadjg_condition_seconds(){ return g_TRACTORHYBRID_saige_getadjg_condition_seconds; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_getadjg_condition_calls(){ return g_TRACTORHYBRID_saige_getadjg_condition_calls; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_getadjg_region_seconds(){ return g_TRACTORHYBRID_saige_getadjg_region_seconds; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_getadjg_region_calls(){ return g_TRACTORHYBRID_saige_getadjg_region_calls; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_getadjg_other_seconds(){ return g_TRACTORHYBRID_saige_getadjg_other_seconds; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_getadjg_other_calls(){ return g_TRACTORHYBRID_saige_getadjg_other_calls; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_spa_seconds(){ return g_TRACTORHYBRID_saige_spa_seconds; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_spa_calls(){ return g_TRACTORHYBRID_saige_spa_calls; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_firth_seconds(){ return g_TRACTORHYBRID_saige_firth_seconds; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_firth_calls(){ return g_TRACTORHYBRID_saige_firth_calls; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_er_seconds(){ return g_TRACTORHYBRID_saige_er_seconds; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_condition_seconds(){ return g_TRACTORHYBRID_saige_condition_seconds; }

// [[Rcpp::export]]
double get_TRACTORHYBRID_saige_region_seconds(){ return g_TRACTORHYBRID_saige_region_seconds; }

namespace SAIGE {
""",
        "get_TRACTORHYBRID_saige_score_seconds",
        "SAIGE_test.cpp internal timing globals",
    )

    text = must_replace(
        text,
        """  	scoreTest(t_GVec, t_Beta, t_seBeta, t_pval_noSPA, pval_noadj, ispvallog, t_altFreq, t_Tstat, t_var1, t_var2, t_gtilde, t_P2Vec, t_gy, is_region, iIndex);
""",
        """  	TRACTORHYBRID_SAIGE_TIMED_SCORE(t_GVec, t_Beta, t_seBeta, t_pval_noSPA, pval_noadj, ispvallog, t_altFreq, t_Tstat, t_var1, t_var2, t_gtilde, t_P2Vec, t_gy, is_region, iIndex);
""",
        "SAIGE_test.cpp scoreTest timer",
        1,
    )
    text = must_replace(
        text,
        """        scoreTestFast(t_GVec, iIndex, t_Beta, t_seBeta, t_pval_noSPA, pval_noadj, ispvallog, t_altFreq, t_Tstat, t_var1, t_var2);
""",
        """        g_TRACTORHYBRID_saige_scorefast_gtilde_target = is_region ? &t_gtilde : nullptr;
        g_TRACTORHYBRID_saige_scorefast_gtilde_ready = false;
        TRACTORHYBRID_SAIGE_TIMED_SCORE_FAST(t_GVec, iIndex, t_Beta, t_seBeta, t_pval_noSPA, pval_noadj, ispvallog, t_altFreq, t_Tstat, t_var1, t_var2);
        if(g_TRACTORHYBRID_saige_scorefast_gtilde_ready){
          is_gtilde = true;
        }
        g_TRACTORHYBRID_saige_scorefast_gtilde_target = nullptr;
""",
        "SAIGE_test.cpp scoreTestFast timer",
        1,
    )
    text = replace_once_if_missing(
        text,
        """    //std::cout << "scoreTestFast " << std::endl;
    arma::vec g1 = t_GVec.elem(t_indexForNonZero);
    arma::mat X1 = m_X.rows(t_indexForNonZero);
    arma::mat A1 = m_XVX_inv_XV.rows(t_indexForNonZero);
    arma::vec mu21;
    arma::vec res1 = m_res.elem(t_indexForNonZero);
    arma::vec Z = A1.t() * g1;
""",
        """    //std::cout << "scoreTestFast " << std::endl;
    double __tractor_saige_scorefast_extract_start = tractor_hybrid_saige_now_seconds();
    arma::vec g1 = t_GVec.elem(t_indexForNonZero);
    arma::mat X1 = m_X.rows(t_indexForNonZero);
    arma::mat A1 = m_XVX_inv_XV.rows(t_indexForNonZero);
    arma::vec mu21;
    arma::vec res1 = m_res.elem(t_indexForNonZero);
    g_TRACTORHYBRID_saige_scorefast_extract_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_scorefast_extract_start;
    double __tractor_saige_scorefast_projection_start = tractor_hybrid_saige_now_seconds();
    arma::vec Z = A1.t() * g1;
""",
        "__tractor_saige_scorefast_extract_start",
        "SAIGE_test.cpp scoreTestFast extract timer",
    )
    text = replace_once_if_missing(
        text,
        """    arma::vec B = X1 * Z;
    arma::vec g1_tilde = g1 - B;
    double var1, var2, S, S1, S2, g1tildemu2;
""",
        """    arma::vec B = X1 * Z;
    arma::vec g1_tilde = g1 - B;
    g_TRACTORHYBRID_saige_scorefast_projection_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_scorefast_projection_start;
    double var1, var2, S, S1, S2, g1tildemu2;
""",
        "g_TRACTORHYBRID_saige_scorefast_projection_seconds +=",
        "SAIGE_test.cpp scoreTestFast projection timer",
    )
    text = replace_once_if_missing(
        text,
        """    //std::cout << "scoreTestFast1 " << std::endl;
    double Bmu2;
    arma::mat  ZtXVXZ = Z.t() * m_XVX * Z;
""",
        """    //std::cout << "scoreTestFast1 " << std::endl;
    double __tractor_saige_scorefast_variance_start = tractor_hybrid_saige_now_seconds();
    double Bmu2;
    arma::mat  ZtXVXZ = Z.t() * m_XVX * Z;
""",
        "__tractor_saige_scorefast_variance_start",
        "SAIGE_test.cpp scoreTestFast variance timer start",
    )
    text = replace_once_if_missing(
        text,
        """    }else if(m_traitType == "quantitative"){
      Bmu2 = dot(g1, B);
      var2 = ZtXVXZ(0,0)*m_tauvec[0] +  dot(g1,g1) - 2*Bmu2;
    }

    var1 = var2 * m_varRatioVal;
""",
        """    }else if(m_traitType == "quantitative"){
      Bmu2 = dot(g1, B);
      var2 = ZtXVXZ(0,0)*m_tauvec[0] +  dot(g1,g1) - 2*Bmu2;
    }
    g_TRACTORHYBRID_saige_scorefast_variance_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_scorefast_variance_start;

    double __tractor_saige_scorefast_result_start = tractor_hybrid_saige_now_seconds();
    var1 = var2 * m_varRatioVal;
""",
        "g_TRACTORHYBRID_saige_scorefast_variance_seconds +=",
        "SAIGE_test.cpp scoreTestFast variance timer end",
    )
    text = replace_once_if_missing(
        text,
        """    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));	
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
}
""",
        """    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));	
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
    g_TRACTORHYBRID_saige_scorefast_result_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_scorefast_result_start;
}
""",
        "g_TRACTORHYBRID_saige_scorefast_result_seconds +=",
        "SAIGE_test.cpp scoreTestFast result timer",
    )
    text = replace_once_if_missing(
        text,
        """  // To increase computational efficiency when lots of GVec elements are 0
 arma::vec m_XVG(m_p, arma::fill::zeros);
  for(int i = 0; i < iIndex.n_elem; i++){
      m_XVG += m_XV.col(iIndex(i)) * t_GVec(iIndex(i));
  }
  g = t_GVec - m_XXVX_inv * m_XVG; 
}
""",
        """  // To increase computational efficiency when lots of GVec elements are 0
  double __tractor_saige_getadjg_accum_start = tractor_hybrid_saige_now_seconds();
 arma::vec m_XVG(m_p, arma::fill::zeros);
  for(int i = 0; i < iIndex.n_elem; i++){
      m_XVG += m_XV.col(iIndex(i)) * t_GVec(iIndex(i));
  }
  g_TRACTORHYBRID_saige_getadjg_accum_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_getadjg_accum_start;
  double __tractor_saige_getadjg_projection_start = tractor_hybrid_saige_now_seconds();
  g = t_GVec - m_XXVX_inv * m_XVG; 
  g_TRACTORHYBRID_saige_getadjg_projection_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_getadjg_projection_start;
}
""",
        "__tractor_saige_getadjg_accum_start",
        "SAIGE_test.cpp getadjGFast internal timers",
    )
    text = replace_once_if_missing(
        text,
        """  unsigned int iIndexComVecSize = iIndexComVec.n_elem;
  unsigned int iIndexSize = iIndex.n_elem; 
  arma::vec gNB(iIndexSize, arma::fill::none);
""",
        """  unsigned int iIndexComVecSize = iIndexComVec.n_elem;
  unsigned int iIndexSize = iIndex.n_elem; 
  double __tractor_saige_alloc_start = tractor_hybrid_saige_now_seconds();
  arma::vec gNB(iIndexSize, arma::fill::none);
""",
        "__tractor_saige_alloc_start",
        "SAIGE_test.cpp allocation timer start",
    )
    text = replace_once_if_missing(
        text,
        """  arma::vec muNA(iIndexComVecSize, arma::fill::none);


  double gmuNB;
""",
        """  arma::vec muNA(iIndexComVecSize, arma::fill::none);
  g_TRACTORHYBRID_saige_alloc_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_alloc_start;


  double gmuNB;
""",
        "g_TRACTORHYBRID_saige_alloc_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_alloc_start;",
        "SAIGE_test.cpp allocation timer end",
    )
    text = text.replace("          getadjGFast(t_GVec, t_gtilde, iIndex);", "          TRACTORHYBRID_SAIGE_TIMED_GETADJG_REASON(g_TRACTORHYBRID_saige_getadjg_spa_seconds, g_TRACTORHYBRID_saige_getadjg_spa_calls, t_GVec, t_gtilde, iIndex);")
    text = text.replace("                getadjGFast(t_GVec, t_gtilde, iIndex);", "                TRACTORHYBRID_SAIGE_TIMED_GETADJG_REASON(g_TRACTORHYBRID_saige_getadjg_firth_seconds, g_TRACTORHYBRID_saige_getadjg_firth_calls, t_GVec, t_gtilde, iIndex);")
    text = text.replace("        \tgetadjGFast(t_GVec, t_gtilde, iIndex);", "        \tTRACTORHYBRID_SAIGE_TIMED_GETADJG_REASON(g_TRACTORHYBRID_saige_getadjg_condition_seconds, g_TRACTORHYBRID_saige_getadjg_condition_calls, t_GVec, t_gtilde, iIndex);")
    text = text.replace("\tgetadjGFast(t_GVec, t_gtilde, iIndex);", "\tTRACTORHYBRID_SAIGE_TIMED_GETADJG_REASON(g_TRACTORHYBRID_saige_getadjg_region_seconds, g_TRACTORHYBRID_saige_getadjg_region_calls, t_GVec, t_gtilde, iIndex);")
    text = replace_once_if_missing(
        text,
        """    getadjGFast(t_GVec, t_gtilde, t_indexForNonZero);
    //getadjG(t_GVec, t_gtilde);


    if(t_is_region && m_traitType == "binary"){
""",
        """    TRACTORHYBRID_SAIGE_TIMED_GETADJG(t_GVec, t_gtilde, t_indexForNonZero);
    //getadjG(t_GVec, t_gtilde);


    if(t_is_region && m_traitType == "binary"){
""",
        "TRACTORHYBRID_SAIGE_TIMED_GETADJG(t_GVec, t_gtilde, t_indexForNonZero);",
        "SAIGE_test.cpp scoreTest getadjGFast timer",
    )
    text = replace_once_if_missing(
        text,
        """void SAIGEClass::extract_anc_stat_for_cond(
        arma::vec & t_GVec,
        arma::vec & t_gtilde,
        arma::mat & t_P2Vec_cond,
        arma::mat & t_VarInvMat_cond,
        arma::mat & t_VarMat_cond,
        arma::vec & t_Tstat_cond,
        arma::uvec & t_indexForNonZero
){
    getadjGFast(t_GVec, t_gtilde, t_indexForNonZero);
    //getadjG(t_GVec, t_gtilde);
""",
        """void SAIGEClass::extract_anc_stat_for_cond(
        arma::vec & t_GVec,
        arma::vec & t_gtilde,
        arma::mat & t_P2Vec_cond,
        arma::mat & t_VarInvMat_cond,
        arma::mat & t_VarMat_cond,
        arma::vec & t_Tstat_cond,
        arma::uvec & t_indexForNonZero
){
    TRACTORHYBRID_SAIGE_TIMED_GETADJG_REASON(g_TRACTORHYBRID_saige_getadjg_condition_seconds, g_TRACTORHYBRID_saige_getadjg_condition_calls, t_GVec, t_gtilde, t_indexForNonZero);
    //getadjG(t_GVec, t_gtilde);
""",
        "TRACTORHYBRID_SAIGE_TIMED_GETADJG_REASON(g_TRACTORHYBRID_saige_getadjg_condition_seconds, g_TRACTORHYBRID_saige_getadjg_condition_calls, t_GVec, t_gtilde, t_indexForNonZero);",
        "SAIGE_test.cpp extract_anc_stat_for_cond getadjGFast timer",
    )
    text = replace_once_if_missing(
        text,
        """	if(!is_gtilde){
                TRACTORHYBRID_SAIGE_TIMED_GETADJG_REASON(g_TRACTORHYBRID_saige_getadjg_spa_seconds, g_TRACTORHYBRID_saige_getadjg_spa_calls, t_GVec, t_gtilde, iIndex);
                is_gtilde = true;
        }
	arma::mat x(t_GVec.n_elem, 2, arma::fill::ones);	
""",
        """	if(!is_gtilde){
                TRACTORHYBRID_SAIGE_TIMED_GETADJG_REASON(g_TRACTORHYBRID_saige_getadjg_firth_seconds, g_TRACTORHYBRID_saige_getadjg_firth_calls, t_GVec, t_gtilde, iIndex);
                is_gtilde = true;
        }
	arma::mat x(t_GVec.n_elem, 2, arma::fill::ones);	
""",
        "TRACTORHYBRID_SAIGE_TIMED_GETADJG_REASON(g_TRACTORHYBRID_saige_getadjg_firth_seconds, g_TRACTORHYBRID_saige_getadjg_firth_calls, t_GVec, t_gtilde, iIndex);",
        "SAIGE_test.cpp Firth getadjGFast reason timer",
    )
    text = text.replace("        \tSPA_fast(m_mu, t_gtilde, q, qinv, pval_noadj, ispvallog, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol1, m_traitType, t_SPApval, t_isSPAConverge);", "        \tTRACTORHYBRID_SAIGE_TIMED_SPA_FAST(m_mu, t_gtilde, q, qinv, pval_noadj, ispvallog, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol1, m_traitType, t_SPApval, t_isSPAConverge);")
    text = text.replace("\t\tSPA(m_mu, t_gtilde, q, qinv, pval_noadj, tol1, ispvallog, m_traitType, t_SPApval, t_isSPAConverge);\t", "\t\tTRACTORHYBRID_SAIGE_TIMED_SPA(m_mu, t_gtilde, q, qinv, pval_noadj, tol1, ispvallog, m_traitType, t_SPApval, t_isSPAConverge);\t")
    text = replace_once_if_missing(
        text,
        """}else{ //if(!t_isER){
double t_GVecMAC = arma::accu(t_GVec); 
""",
        """}else{ //if(!t_isER){
double __tractor_saige_er_start = tractor_hybrid_saige_now_seconds();
double t_GVecMAC = arma::accu(t_GVec); 
""",
        "__tractor_saige_er_start",
        "SAIGE_test.cpp ER timer start",
    )
    text = replace_once_if_missing(
        text,
        """    }\t    
}

   if(t_isFirth){
""",
        """    }\t    
  g_TRACTORHYBRID_saige_er_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_er_start;
}

   if(t_isFirth){
   double __tractor_saige_firth_start = tractor_hybrid_saige_now_seconds();
   g_TRACTORHYBRID_saige_firth_calls += 1.0;
""",
        "__tractor_saige_firth_start",
        "SAIGE_test.cpp ER timer end and Firth timer start",
    )
    text = replace_once_if_missing(
        text,
        """\t//back calculates se based on beta from firth adjustion and the p-value that accounts for case-control imbalance
\tt_seBeta = fabs(t_Beta)/fabs(t_qval_Firth);
   }
   
 //arma::vec timeoutput4 = getTime();
""",
        """\t//back calculates se based on beta from firth adjustion and the p-value that accounts for case-control imbalance
\tt_seBeta = fabs(t_Beta)/fabs(t_qval_Firth);
   g_TRACTORHYBRID_saige_firth_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_firth_start;
   }
   
 //arma::vec timeoutput4 = getTime();
""",
        "g_TRACTORHYBRID_saige_firth_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_firth_start;",
        "SAIGE_test.cpp Firth timer end",
    )
    text = replace_once_if_missing(
        text,
        """   if(t_isCondition){
\tif(!is_gtilde){
""",
        """   if(t_isCondition){
   double __tractor_saige_condition_start = tractor_hybrid_saige_now_seconds();
\tif(!is_gtilde){
""",
        "__tractor_saige_condition_start",
        "SAIGE_test.cpp condition timer start",
    )
    text = replace_once_if_missing(
        text,
        """ }

    gNA.clear();
""",
        """   g_TRACTORHYBRID_saige_condition_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_condition_start;
 }

    gNA.clear();
""",
        "g_TRACTORHYBRID_saige_condition_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_condition_start;",
        "SAIGE_test.cpp condition timer end",
    )
    text = replace_once_if_missing(
        text,
        """    if(is_region && isScoreFast){

      t_gy = dot(t_gtilde, m_y);
""",
        """    if(is_region && isScoreFast){
      double __tractor_saige_region_start = tractor_hybrid_saige_now_seconds();

      t_gy = dot(t_gtilde, m_y);
""",
        "__tractor_saige_region_start",
        "SAIGE_test.cpp region timer start",
    )
    text = replace_once_if_missing(
        text,
        """      }
    }
}
""",
        """      }
      g_TRACTORHYBRID_saige_region_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_region_start;
    }
}
""",
        "g_TRACTORHYBRID_saige_region_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_region_start;",
        "SAIGE_test.cpp region timer end",
    )
    text = replace_between(
        text,
        "void SAIGEClass::scoreTestFast(",
        "\n\n\nvoid SAIGEClass::getadjG(",
        """void SAIGEClass::scoreTestFast(arma::vec & t_GVec,
                     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
\t\t     double& t_pval,
                     bool& t_islogp,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2){
    // TRACTORHYBRID optimized sparse score path: same algebra as the original
    // scoreTestFast, but avoids large row-subset temporaries and scans
    // Armadillo's column-major matrices by column.
    const arma::uword p = m_X.n_cols;

    double __tractor_saige_scorefast_extract_start = tractor_hybrid_saige_now_seconds();
    arma::vec Z(p, arma::fill::zeros);
    for(arma::uword k = 0; k < p; ++k){
      const double* avec = m_XVX_inv_XV.colptr(k);
      double z = 0.0;
      for(arma::uword ii = 0; ii < t_indexForNonZero.n_elem; ++ii){
        const arma::uword row = t_indexForNonZero(ii);
        z += avec[row] * t_GVec(row);
      }
      Z(k) = z;
    }
    g_TRACTORHYBRID_saige_scorefast_extract_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_scorefast_extract_start;

    double __tractor_saige_scorefast_gtilde_start = tractor_hybrid_saige_now_seconds();
    arma::vec* __tractor_scorefast_gtilde = g_TRACTORHYBRID_saige_scorefast_gtilde_target;
    double* __tractor_scorefast_gtilde_ptr = nullptr;
    if(__tractor_scorefast_gtilde != nullptr){
      __tractor_scorefast_gtilde->set_size(t_GVec.n_elem);
      __tractor_scorefast_gtilde->zeros();
      __tractor_scorefast_gtilde_ptr = __tractor_scorefast_gtilde->memptr();
      const arma::uword n = t_GVec.n_elem;
      for(arma::uword k = 0; k < p; ++k){
        const double zk = Z(k);
        if(zk == 0.0){
          continue;
        }
        const double* xvec = m_X.colptr(k);
        for(arma::uword row = 0; row < n; ++row){
          __tractor_scorefast_gtilde_ptr[row] -= xvec[row] * zk;
        }
      }
      for(arma::uword ii = 0; ii < t_indexForNonZero.n_elem; ++ii){
        const arma::uword row = t_indexForNonZero(ii);
        __tractor_scorefast_gtilde_ptr[row] += t_GVec(row);
      }
      g_TRACTORHYBRID_saige_scorefast_gtilde_ready = true;
      g_TRACTORHYBRID_saige_scorefast_gtilde_calls += 1.0;
      g_TRACTORHYBRID_saige_scorefast_gtilde_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_scorefast_gtilde_start;
    }

    double __tractor_saige_scorefast_projection_start = tractor_hybrid_saige_now_seconds();
    arma::vec res1X1(p, arma::fill::zeros);
    double S1 = 0.0;
    double Bmu2 = 0.0;
    double g1tildemu2 = 0.0;
    double g1g1 = 0.0;

    if(__tractor_scorefast_gtilde_ptr != nullptr){
      for(arma::uword k = 0; k < p; ++k){
        const double* xvec = m_X.colptr(k);
        double rx = 0.0;
        for(arma::uword ii = 0; ii < t_indexForNonZero.n_elem; ++ii){
          const arma::uword row = t_indexForNonZero(ii);
          rx += m_res(row) * xvec[row];
        }
        res1X1(k) = rx;
      }

      for(arma::uword ii = 0; ii < t_indexForNonZero.n_elem; ++ii){
        const arma::uword row = t_indexForNonZero(ii);
        const double gv = t_GVec(row);
        const double res = m_res(row);
        const double g1_tilde = __tractor_scorefast_gtilde_ptr[row];
        const double B = gv - g1_tilde;

        S1 += res * g1_tilde;

        if(m_traitType == "binary" || m_traitType == "survival"){
          const double mu2 = m_mu2(row);
          g1tildemu2 += g1_tilde * g1_tilde * mu2;
          Bmu2 += B * B * mu2;
        }else if(m_traitType == "quantitative"){
          Bmu2 += gv * B;
          g1g1 += gv * gv;
        }
      }
    }else{
      arma::vec Bvec(t_indexForNonZero.n_elem, arma::fill::zeros);
      for(arma::uword k = 0; k < p; ++k){
        const double* xvec = m_X.colptr(k);
        const double zk = Z(k);
        double rx = 0.0;
        for(arma::uword ii = 0; ii < t_indexForNonZero.n_elem; ++ii){
          const arma::uword row = t_indexForNonZero(ii);
          const double x = xvec[row];
          Bvec(ii) += x * zk;
          rx += m_res(row) * x;
        }
        res1X1(k) = rx;
      }

      for(arma::uword ii = 0; ii < t_indexForNonZero.n_elem; ++ii){
        const arma::uword row = t_indexForNonZero(ii);
        const double gv = t_GVec(row);
        const double res = m_res(row);
        const double B = Bvec(ii);

        const double g1_tilde = gv - B;
        S1 += res * g1_tilde;

        if(m_traitType == "binary" || m_traitType == "survival"){
          const double mu2 = m_mu2(row);
          g1tildemu2 += g1_tilde * g1_tilde * mu2;
          Bmu2 += B * B * mu2;
        }else if(m_traitType == "quantitative"){
          Bmu2 += gv * B;
          g1g1 += gv * gv;
        }
      }
    }
    g_TRACTORHYBRID_saige_scorefast_projection_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_scorefast_projection_start;

    double __tractor_saige_scorefast_variance_start = tractor_hybrid_saige_now_seconds();
    const double ZtXVXZ = arma::as_scalar(Z.t() * m_XVX * Z);
    double var2 = 0.0;
    if(m_traitType == "binary" || m_traitType == "survival"){
      var2 = ZtXVXZ - Bmu2 + g1tildemu2;
    }else if(m_traitType == "quantitative"){
      var2 = ZtXVXZ * m_tauvec[0] + g1g1 - 2 * Bmu2;
    }
    g_TRACTORHYBRID_saige_scorefast_variance_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_scorefast_variance_start;

    double __tractor_saige_scorefast_result_start = tractor_hybrid_saige_now_seconds();
    const double var1 = var2 * m_varRatioVal;
    double S2 = 0.0;
    for(arma::uword k = 0; k < p; ++k){
      S2 -= (m_S_a(k) - res1X1(k)) * Z(k);
    }
    const double S = (S1 + S2) / m_tauvec[0];

    double stat = S*S/var1;
    if (var1 <= std::numeric_limits<double>::min()){
          t_pval = 1;
    }else{
      if(!std::isnan(stat) && std::isfinite(stat)){
          boost::math::chi_squared chisq_dist(1);
          t_pval = boost::math::cdf(complement(chisq_dist, stat));

      }else{
          t_pval = 1;
\t  stat = 0.0;
      }
    }
    char pValueBuf[100];
    if (t_pval != 0){
        sprintf(pValueBuf, "%.6E", t_pval);
\tt_islogp = false;
    }else {
\tdouble logp = R::pchisq(stat,1,false,true);
        double log10p = logp/(log(10));
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
         }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
\tt_pval = logp;
\tt_islogp = true;
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr;
    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
    g_TRACTORHYBRID_saige_scorefast_result_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_scorefast_result_start;
}
""",
        "SAIGE_test.cpp optimized scoreTestFast replacement",
    )
    text = replace_between(
        text,
        "void SAIGEClass::getadjGFast(",
        "\n\n\nvoid SAIGEClass::get_mu(",
        """void SAIGEClass::getadjGFast(arma::vec & t_GVec, arma::vec & g, arma::uvec & iIndex)
{

  // To increase computational efficiency when lots of GVec elements are 0
  double __tractor_saige_getadjg_accum_start = tractor_hybrid_saige_now_seconds();
  arma::vec m_XVG(m_p, arma::fill::zeros);
  for(arma::uword i = 0; i < iIndex.n_elem; ++i){
      const arma::uword col = iIndex(i);
      const double gv = t_GVec(col);
      for(arma::uword k = 0; k < m_XVG.n_elem; ++k){
          m_XVG(k) += m_XV(k, col) * gv;
      }
  }
  g_TRACTORHYBRID_saige_getadjg_accum_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_getadjg_accum_start;

  double __tractor_saige_getadjg_projection_start = tractor_hybrid_saige_now_seconds();
  g = t_GVec;
  double* gptr = g.memptr();
  const arma::uword n = g.n_elem;
  for(arma::uword k = 0; k < m_XVG.n_elem; ++k){
      const double coef = m_XVG(k);
      if(coef == 0.0){
          continue;
      }
      const double* xptr = m_XXVX_inv.colptr(k);
      for(arma::uword row = 0; row < n; ++row){
          gptr[row] -= xptr[row] * coef;
      }
  }
  g_TRACTORHYBRID_saige_getadjg_projection_seconds += tractor_hybrid_saige_now_seconds() - __tractor_saige_getadjg_projection_start;
}
""",
        "SAIGE_test.cpp optimized getadjGFast replacement",
    )
    write(path, text)


def patch_geno_r():
    path = APP / "R" / "Geno.R"
    text = read(path)
    text = text.replace("sampleInModel = NULL){", "sampleInModel = NULL,\n                 tractorHybridPrefix = \"\"){", 1)
    text = text.replace("sampleInModel = NULL)\n{", "sampleInModel = NULL,\n\t\t tractorHybridPrefix = \"\")\n{", 1)

    text = replace_once_if_missing(
        text,
        "\n    }else if(bedFile != \"\"){\n",
        """
    }else if(tractorHybridPrefix != ""){
        Check_File_Exist(paste0(tractorHybridPrefix, ".meta"), "tractorHybridPrefix meta")
        Check_File_Exist(paste0(tractorHybridPrefix, ".samples"), "tractorHybridPrefix samples")
        Check_File_Exist(paste0(tractorHybridPrefix, ".common.variant.mks"), "tractorHybridPrefix common mks")
        Check_File_Exist(paste0(tractorHybridPrefix, ".rare.variant.mks"), "tractorHybridPrefix rare mks")
        Check_File_Exist(paste0(tractorHybridPrefix, ".common.geno.bin"), "tractorHybridPrefix common bin")
        Check_File_Exist(paste0(tractorHybridPrefix, ".rare.carrier.bin"), "tractorHybridPrefix rare bin")
        Check_File_Exist(paste0(tractorHybridPrefix, ".ancblock.bin"), "tractorHybridPrefix anc bin")
        Check_File_Exist(paste0(tractorHybridPrefix, ".ancblock.mks"), "tractorHybridPrefix anc mks")
        Check_File_Exist(paste0(tractorHybridPrefix, ".ancblock.idx"), "tractorHybridPrefix anc idx")
        dosageFileType = "tractor_hybrid"
    }else if(bedFile != ""){
""",
        'dosageFileType = "tractor_hybrid"',
        "Geno.R checkGenoInput branch",
    )

    text = text.replace(
        "famFile = famFile,\n\t\t sampleInModel = sampleInModel)",
        "famFile = famFile,\n\t\t sampleInModel = sampleInModel,\n                 tractorHybridPrefix = tractorHybridPrefix)",
        1,
    )

    text = replace_once_if_missing(
        text,
        """  if(dosageFileType == "vcf"){
    if(idstoIncludeFile != "" & rangestoIncludeFile != ""){
      stop("We currently do not support both 'idstoIncludeFile' and 'rangestoIncludeFile' at the same time for vcf files\\n")
    }
    #if(chrom==""){
    #  stop("chrom needs to be specified for VCF/BCF/SAV input\\n")
    #}
    markerInfo = NULL
  }
""",
        """  if(dosageFileType == "vcf"){
    if(idstoIncludeFile != "" & rangestoIncludeFile != ""){
      stop("We currently do not support both 'idstoIncludeFile' and 'rangestoIncludeFile' at the same time for vcf files\\n")
    }
    #if(chrom==""){
    #  stop("chrom needs to be specified for VCF/BCF/SAV input\\n")
    #}
    markerInfo = NULL
  }

  if(dosageFileType == "tractor_hybrid"){
    if(idstoIncludeFile != "" | rangestoIncludeFile != ""){
      stop("tractor_hybrid streaming step2 currently supports chromosome traversal only; do not set idstoIncludeFile or rangestoIncludeFile\\n")
    }
    markerInfo = NULL
    setTRACTORHYBRIDobjInCPP(tractorHybridPrefix, sampleInModel)
  }
""",
        "setTRACTORHYBRIDobjInCPP",
        "Geno.R setGenoInput branch",
    )

    write(path, text)


def patch_saige_test_main():
    path = APP / "R" / "SAIGE_Test_main.R"
    text = read(path)
    text = text.replace('famFile="",\n                 AlleleOrder = "alt-first"', 'famFile="",\n                 tractorHybridPrefix = "",\n                 AlleleOrder = "alt-first"', 1)
    text = text.replace("famFile=famFile,\n                 idstoIncludeFile", "famFile=famFile,\n                 tractorHybridPrefix = tractorHybridPrefix,\n                 idstoIncludeFile", 1)
    text = replace_once_if_missing(
        text,
        """    if (isCondition) {
        n = length(obj.model$y) #sample size
""",
        """    if (isCondition) {
        if(genoType == "tractor_hybrid"){
            stop("tractor_hybrid streaming step2 currently does not support conditional analysis by marker ID")
        }
        n = length(obj.model$y) #sample size
""",
        'tractor_hybrid streaming step2 currently does not support conditional analysis',
        "SAIGE_Test_main.R condition guard",
    )
    write(path, text)


def patch_spatest_tractor():
    path = APP / "R" / "SAIGE_SPATest_Tractor.R"
    text = read(path)
    text = text.replace('if(genoType != "vcf"){', 'if(genoType != "vcf" && genoType != "tractor_hybrid"){', 1)
    old = """  }else{
   if(!isAnyInclude){
    if(chrom == ""){
      stop("chrom needs to be specified for single-variant assoc tests when using VCF as input\\n")
    }else{
      set_iterator_inVcf("", chrom, 1, 250000000)
    }
   }
    if(outIndex > 1){
        move_forward_iterator_Vcf(outIndex*nMarkersEachChunk)
    }
    isVcfEnd =  check_Vcf_end()
    if(!isVcfEnd){
        #outIndex = 1
        genoIndex = rep("0", nMarkersEachChunk)
        genoIndex_prev = rep("0", nMarkersEachChunk)
        #nChunks = outIndex + 1
        is_marker_test = TRUE
        i = outIndex
    }else{
        is_marker_test = FALSE
        stop("No markers are left in VCF")
    }
  }
"""
    new = """  }else{
   if(!isAnyInclude){
    if(chrom == ""){
      stop("chrom needs to be specified for single-variant assoc tests when using VCF or tractor_hybrid as input\\n")
    }else if(genoType == "vcf"){
      set_iterator_inVcf("", chrom, 1, 250000000)
    }else if(genoType == "tractor_hybrid"){
      set_iterator_inTRACTORHYBRID(chrom, 1, 250000000)
    }
   }
    if(outIndex > 1){
      if(genoType == "vcf"){
        move_forward_iterator_Vcf(outIndex*nMarkersEachChunk)
      }else if(genoType == "tractor_hybrid"){
        move_forward_iterator_TRACTORHYBRID(outIndex*nMarkersEachChunk)
      }
    }
    isStreamEnd = if(genoType == "vcf") check_Vcf_end() else check_TRACTORHYBRID_end()
    if(!isStreamEnd){
        #outIndex = 1
        genoIndex = rep("0", nMarkersEachChunk)
        genoIndex_prev = rep("0", nMarkersEachChunk)
        #nChunks = outIndex + 1
        is_marker_test = TRUE
        i = outIndex
    }else{
        is_marker_test = FALSE
        stop("No markers are left in genotype stream")
    }
  }
"""
    text = must_replace(text, old, new, "SAIGE_SPATest_Tractor.R streaming setup", 1)
    text = text.replace('if(genoType != "vcf"){\n      tempList = genoIndexList[[i]]', 'if(genoType != "vcf" && genoType != "tractor_hybrid"){\n      tempList = genoIndexList[[i]]', 1)
    text = text.replace(
        'if(genoType != "vcf"){\n      cat(paste0("(",Sys.time(),") ---- Analyzing Chunk ", i, "/", nChunks, ": chrom ", chrom," ---- \\\\n"))',
        'if(genoType != "vcf" && genoType != "tractor_hybrid"){\n      cat(paste0("(",Sys.time(),") ---- Analyzing Chunk ", i, "/", nChunks, ": chrom ", chrom," ---- \\\\n"))',
        1,
    )
    text = text.replace(
        'if(genoType != "vcf"){\n      cat(paste0("(",Sys.time(),") ---- Analyzing Chunk ',
        'if(genoType != "vcf" && genoType != "tractor_hybrid"){\n      cat(paste0("(",Sys.time(),") ---- Analyzing Chunk ',
        1,
    )
    text = must_replace(
        text,
        "  mainMarkerAdmixedInCPP(genoType, traitType, genoIndex_prev, genoIndex, isMoreOutput, isImputation, isFirth, NumberofANC, pvalcutoff_of_haplotype)\n",
        """  if(genoType == "tractor_hybrid"){
    .tractor_io_before = get_TRACTORHYBRID_io_seconds()
    .tractor_decode_before = get_TRACTORHYBRID_decode_seconds()
    .tractor_get_marker_before = get_TRACTORHYBRID_get_marker_seconds()
    .tractor_pvalue_before = get_TRACTORHYBRID_pvalue_seconds()
    .tractor_output_before = get_TRACTORHYBRID_output_seconds()
    .tractor_reset_before = get_TRACTORHYBRID_reset_seconds()
    .tractor_condition_before = get_TRACTORHYBRID_condition_seconds()
    .tractor_impute_qc_before = get_TRACTORHYBRID_impute_qc_seconds()
    .tractor_joint_before = get_TRACTORHYBRID_joint_seconds()
    .tractor_variance_before = get_TRACTORHYBRID_variance_seconds()
    .tractor_condition_cache_hits_before = get_TRACTORHYBRID_condition_cache_hits()
    .tractor_condition_cache_misses_before = get_TRACTORHYBRID_condition_cache_misses()
    .tractor_saige_score_before = get_TRACTORHYBRID_saige_score_seconds()
    .tractor_saige_score_fast_calls_before = get_TRACTORHYBRID_saige_score_fast_calls()
    .tractor_saige_score_slow_calls_before = get_TRACTORHYBRID_saige_score_slow_calls()
    .tractor_saige_scorefast_extract_before = get_TRACTORHYBRID_saige_scorefast_extract_seconds()
    .tractor_saige_scorefast_projection_before = get_TRACTORHYBRID_saige_scorefast_projection_seconds()
    .tractor_saige_scorefast_variance_before = get_TRACTORHYBRID_saige_scorefast_variance_seconds()
    .tractor_saige_scorefast_result_before = get_TRACTORHYBRID_saige_scorefast_result_seconds()
    .tractor_saige_scorefast_gtilde_before = get_TRACTORHYBRID_saige_scorefast_gtilde_seconds()
    .tractor_saige_scorefast_gtilde_calls_before = get_TRACTORHYBRID_saige_scorefast_gtilde_calls()
    .tractor_saige_alloc_before = get_TRACTORHYBRID_saige_alloc_seconds()
    .tractor_saige_getadjg_before = get_TRACTORHYBRID_saige_getadjg_seconds()
    .tractor_saige_getadjg_calls_before = get_TRACTORHYBRID_saige_getadjg_calls()
    .tractor_saige_getadjg_accum_before = get_TRACTORHYBRID_saige_getadjg_accum_seconds()
    .tractor_saige_getadjg_projection_before = get_TRACTORHYBRID_saige_getadjg_projection_seconds()
    .tractor_saige_getadjg_spa_before = get_TRACTORHYBRID_saige_getadjg_spa_seconds()
    .tractor_saige_getadjg_spa_calls_before = get_TRACTORHYBRID_saige_getadjg_spa_calls()
    .tractor_saige_getadjg_firth_before = get_TRACTORHYBRID_saige_getadjg_firth_seconds()
    .tractor_saige_getadjg_firth_calls_before = get_TRACTORHYBRID_saige_getadjg_firth_calls()
    .tractor_saige_getadjg_condition_before = get_TRACTORHYBRID_saige_getadjg_condition_seconds()
    .tractor_saige_getadjg_condition_calls_before = get_TRACTORHYBRID_saige_getadjg_condition_calls()
    .tractor_saige_getadjg_region_before = get_TRACTORHYBRID_saige_getadjg_region_seconds()
    .tractor_saige_getadjg_region_calls_before = get_TRACTORHYBRID_saige_getadjg_region_calls()
    .tractor_saige_getadjg_other_before = get_TRACTORHYBRID_saige_getadjg_other_seconds()
    .tractor_saige_getadjg_other_calls_before = get_TRACTORHYBRID_saige_getadjg_other_calls()
    .tractor_saige_spa_before = get_TRACTORHYBRID_saige_spa_seconds()
    .tractor_saige_spa_calls_before = get_TRACTORHYBRID_saige_spa_calls()
    .tractor_saige_firth_before = get_TRACTORHYBRID_saige_firth_seconds()
    .tractor_saige_firth_calls_before = get_TRACTORHYBRID_saige_firth_calls()
    .tractor_saige_er_before = get_TRACTORHYBRID_saige_er_seconds()
    .tractor_saige_condition_before = get_TRACTORHYBRID_saige_condition_seconds()
    .tractor_saige_region_before = get_TRACTORHYBRID_saige_region_seconds()
    .tractor_chunk_time = system.time({
      mainMarkerAdmixedInCPP(genoType, traitType, genoIndex_prev, genoIndex, isMoreOutput, isImputation, isFirth, NumberofANC, pvalcutoff_of_haplotype)
    })
    .tractor_io_delta = get_TRACTORHYBRID_io_seconds() - .tractor_io_before
    .tractor_decode_delta = get_TRACTORHYBRID_decode_seconds() - .tractor_decode_before
    .tractor_get_marker_delta = get_TRACTORHYBRID_get_marker_seconds() - .tractor_get_marker_before
    .tractor_pvalue_delta = get_TRACTORHYBRID_pvalue_seconds() - .tractor_pvalue_before
    .tractor_output_delta = get_TRACTORHYBRID_output_seconds() - .tractor_output_before
    .tractor_reset_delta = get_TRACTORHYBRID_reset_seconds() - .tractor_reset_before
    .tractor_condition_delta = get_TRACTORHYBRID_condition_seconds() - .tractor_condition_before
    .tractor_impute_qc_delta = get_TRACTORHYBRID_impute_qc_seconds() - .tractor_impute_qc_before
    .tractor_joint_delta = get_TRACTORHYBRID_joint_seconds() - .tractor_joint_before
    .tractor_variance_delta = get_TRACTORHYBRID_variance_seconds() - .tractor_variance_before
    .tractor_condition_cache_hits_delta = get_TRACTORHYBRID_condition_cache_hits() - .tractor_condition_cache_hits_before
    .tractor_condition_cache_misses_delta = get_TRACTORHYBRID_condition_cache_misses() - .tractor_condition_cache_misses_before
    .tractor_saige_score_delta = get_TRACTORHYBRID_saige_score_seconds() - .tractor_saige_score_before
    .tractor_saige_score_fast_calls_delta = get_TRACTORHYBRID_saige_score_fast_calls() - .tractor_saige_score_fast_calls_before
    .tractor_saige_score_slow_calls_delta = get_TRACTORHYBRID_saige_score_slow_calls() - .tractor_saige_score_slow_calls_before
    .tractor_saige_scorefast_extract_delta = get_TRACTORHYBRID_saige_scorefast_extract_seconds() - .tractor_saige_scorefast_extract_before
    .tractor_saige_scorefast_projection_delta = get_TRACTORHYBRID_saige_scorefast_projection_seconds() - .tractor_saige_scorefast_projection_before
    .tractor_saige_scorefast_variance_delta = get_TRACTORHYBRID_saige_scorefast_variance_seconds() - .tractor_saige_scorefast_variance_before
    .tractor_saige_scorefast_result_delta = get_TRACTORHYBRID_saige_scorefast_result_seconds() - .tractor_saige_scorefast_result_before
    .tractor_saige_scorefast_gtilde_delta = get_TRACTORHYBRID_saige_scorefast_gtilde_seconds() - .tractor_saige_scorefast_gtilde_before
    .tractor_saige_scorefast_gtilde_calls_delta = get_TRACTORHYBRID_saige_scorefast_gtilde_calls() - .tractor_saige_scorefast_gtilde_calls_before
    .tractor_saige_alloc_delta = get_TRACTORHYBRID_saige_alloc_seconds() - .tractor_saige_alloc_before
    .tractor_saige_getadjg_delta = get_TRACTORHYBRID_saige_getadjg_seconds() - .tractor_saige_getadjg_before
    .tractor_saige_getadjg_calls_delta = get_TRACTORHYBRID_saige_getadjg_calls() - .tractor_saige_getadjg_calls_before
    .tractor_saige_getadjg_accum_delta = get_TRACTORHYBRID_saige_getadjg_accum_seconds() - .tractor_saige_getadjg_accum_before
    .tractor_saige_getadjg_projection_delta = get_TRACTORHYBRID_saige_getadjg_projection_seconds() - .tractor_saige_getadjg_projection_before
    .tractor_saige_getadjg_spa_delta = get_TRACTORHYBRID_saige_getadjg_spa_seconds() - .tractor_saige_getadjg_spa_before
    .tractor_saige_getadjg_spa_calls_delta = get_TRACTORHYBRID_saige_getadjg_spa_calls() - .tractor_saige_getadjg_spa_calls_before
    .tractor_saige_getadjg_firth_delta = get_TRACTORHYBRID_saige_getadjg_firth_seconds() - .tractor_saige_getadjg_firth_before
    .tractor_saige_getadjg_firth_calls_delta = get_TRACTORHYBRID_saige_getadjg_firth_calls() - .tractor_saige_getadjg_firth_calls_before
    .tractor_saige_getadjg_condition_delta = get_TRACTORHYBRID_saige_getadjg_condition_seconds() - .tractor_saige_getadjg_condition_before
    .tractor_saige_getadjg_condition_calls_delta = get_TRACTORHYBRID_saige_getadjg_condition_calls() - .tractor_saige_getadjg_condition_calls_before
    .tractor_saige_getadjg_region_delta = get_TRACTORHYBRID_saige_getadjg_region_seconds() - .tractor_saige_getadjg_region_before
    .tractor_saige_getadjg_region_calls_delta = get_TRACTORHYBRID_saige_getadjg_region_calls() - .tractor_saige_getadjg_region_calls_before
    .tractor_saige_getadjg_other_delta = get_TRACTORHYBRID_saige_getadjg_other_seconds() - .tractor_saige_getadjg_other_before
    .tractor_saige_getadjg_other_calls_delta = get_TRACTORHYBRID_saige_getadjg_other_calls() - .tractor_saige_getadjg_other_calls_before
    .tractor_saige_spa_delta = get_TRACTORHYBRID_saige_spa_seconds() - .tractor_saige_spa_before
    .tractor_saige_spa_calls_delta = get_TRACTORHYBRID_saige_spa_calls() - .tractor_saige_spa_calls_before
    .tractor_saige_firth_delta = get_TRACTORHYBRID_saige_firth_seconds() - .tractor_saige_firth_before
    .tractor_saige_firth_calls_delta = get_TRACTORHYBRID_saige_firth_calls() - .tractor_saige_firth_calls_before
    .tractor_saige_er_delta = get_TRACTORHYBRID_saige_er_seconds() - .tractor_saige_er_before
    .tractor_saige_condition_delta = get_TRACTORHYBRID_saige_condition_seconds() - .tractor_saige_condition_before
    .tractor_saige_region_delta = get_TRACTORHYBRID_saige_region_seconds() - .tractor_saige_region_before
    .tractor_total = as.numeric(.tractor_chunk_time[["elapsed"]])
    .tractor_get_marker_overhead = max(0, .tractor_get_marker_delta - .tractor_io_delta - .tractor_decode_delta)
    .tractor_other = max(0, .tractor_total - .tractor_get_marker_delta - .tractor_pvalue_delta - .tractor_output_delta - .tractor_reset_delta - .tractor_impute_qc_delta - .tractor_joint_delta - .tractor_variance_delta)
    cat(sprintf("tractor_hybrid timing chunk %s: total=%.3fs input_io=%.3fs reader_decode=%.3fs hybrid_get_marker=%.3fs get_marker_overhead=%.3fs marker_pvalue=%.3fs reset_zero=%.3fs condition_total=%.3fs condition_cache_hit=%.0f condition_cache_miss=%.0f impute_qc=%.3fs variance_ratio=%.3fs joint_cct_spa=%.3fs output_write=%.3fs other_cpp_or_r=%.3fs\\n",
                i, .tractor_total, .tractor_io_delta, .tractor_decode_delta, .tractor_get_marker_delta, .tractor_get_marker_overhead, .tractor_pvalue_delta, .tractor_reset_delta, .tractor_condition_delta, .tractor_condition_cache_hits_delta, .tractor_condition_cache_misses_delta, .tractor_impute_qc_delta, .tractor_variance_delta, .tractor_joint_delta, .tractor_output_delta, .tractor_other))
    cat(sprintf("tractor_hybrid pvalue timing chunk %s: pvalue_score=%.3fs score_fast_calls=%.0f score_slow_calls=%.0f scorefast_extract=%.3fs scorefast_projection=%.3fs scorefast_variance=%.3fs scorefast_result=%.3fs scorefast_gtilde=%.3fs scorefast_gtilde_calls=%.0f pvalue_alloc=%.3fs pvalue_getadjg=%.3fs getadjg_calls=%.0f getadjg_accum=%.3fs getadjg_projection=%.3fs getadjg_spa=%.3fs getadjg_spa_calls=%.0f getadjg_firth=%.3fs getadjg_firth_calls=%.0f getadjg_condition=%.3fs getadjg_condition_calls=%.0f getadjg_region=%.3fs getadjg_region_calls=%.0f getadjg_other=%.3fs getadjg_other_calls=%.0f pvalue_spa=%.3fs spa_calls=%.0f pvalue_firth=%.3fs firth_calls=%.0f pvalue_er=%.3fs pvalue_condition_adjust=%.3fs pvalue_region_finalize=%.3fs\\n",
                i, .tractor_saige_score_delta, .tractor_saige_score_fast_calls_delta, .tractor_saige_score_slow_calls_delta, .tractor_saige_scorefast_extract_delta, .tractor_saige_scorefast_projection_delta, .tractor_saige_scorefast_variance_delta, .tractor_saige_scorefast_result_delta, .tractor_saige_scorefast_gtilde_delta, .tractor_saige_scorefast_gtilde_calls_delta, .tractor_saige_alloc_delta, .tractor_saige_getadjg_delta, .tractor_saige_getadjg_calls_delta, .tractor_saige_getadjg_accum_delta, .tractor_saige_getadjg_projection_delta, .tractor_saige_getadjg_spa_delta, .tractor_saige_getadjg_spa_calls_delta, .tractor_saige_getadjg_firth_delta, .tractor_saige_getadjg_firth_calls_delta, .tractor_saige_getadjg_condition_delta, .tractor_saige_getadjg_condition_calls_delta, .tractor_saige_getadjg_region_delta, .tractor_saige_getadjg_region_calls_delta, .tractor_saige_getadjg_other_delta, .tractor_saige_getadjg_other_calls_delta, .tractor_saige_spa_delta, .tractor_saige_spa_calls_delta, .tractor_saige_firth_delta, .tractor_saige_firth_calls_delta, .tractor_saige_er_delta, .tractor_saige_condition_delta, .tractor_saige_region_delta))
  }else{
    mainMarkerAdmixedInCPP(genoType, traitType, genoIndex_prev, genoIndex, isMoreOutput, isImputation, isFirth, NumberofANC, pvalcutoff_of_haplotype)
  }
""",
        "SAIGE_SPATest_Tractor.R timing wrapper",
        1,
    )
    text = must_replace(
        text,
        """  if(genoType == "vcf"){
    isEnd_Output =  check_Vcf_end()
  }else{
    isEnd_Output = (i==nChunks)
  }
""",
        """  if(genoType == "vcf"){
    isEnd_Output =  check_Vcf_end()
  }else if(genoType == "tractor_hybrid"){
    isEnd_Output = check_TRACTORHYBRID_end()
  }else{
    isEnd_Output = (i==nChunks)
  }
""",
        "SAIGE_SPATest_Tractor.R end output",
        1,
    )
    text = must_replace(
        text,
        """  if(genoType == "vcf"){
    isVcfEnd =  check_Vcf_end()
    cat("isVcfEnd ", isVcfEnd, "\\n")
    if(isVcfEnd){
        is_marker_test = FALSE
    }
  }else{
    if(i > nChunks){
      is_marker_test = FALSE
    }
  }
""",
        """  if(genoType == "vcf"){
    isStreamEnd =  check_Vcf_end()
    cat("isVcfEnd ", isStreamEnd, "\\n")
    if(isStreamEnd){
        is_marker_test = FALSE
    }
  }else if(genoType == "tractor_hybrid"){
    isStreamEnd = check_TRACTORHYBRID_end()
    cat("isTRACTORHYBRIDEnd ", isStreamEnd, "\\n")
    if(isStreamEnd){
        is_marker_test = FALSE
    }
  }else{
    if(i > nChunks){
      is_marker_test = FALSE
    }
  }
""",
        "SAIGE_SPATest_Tractor.R continuation",
        1,
    )
    write(path, text)


def patch_spatest_marker():
    path = APP / "R" / "SAIGE_SPATest_Marker.R"
    text = read(path)
    insert_after_vcf = """  }else if(genoType == "bgen"){
"""
    hybrid_branch = """  }else if(genoType == "tractor_hybrid"){
   if(!isAnyInclude){
    if(chrom == ""){
      stop("chrom needs to be specified for single-variant assoc tests when using tractor_hybrid as input\\n")
    }else{
      set_iterator_inTRACTORHYBRID(chrom, 1, 250000000)
    }
   }
    if(outIndex > 1){
      move_forward_iterator_TRACTORHYBRID(outIndex*nMarkersEachChunk)
    }
    isHybridEnd = check_TRACTORHYBRID_end()
    if(!isHybridEnd){
      genoIndex = rep("0", nMarkersEachChunk)
      genoIndex_prev = rep("0", nMarkersEachChunk)
      is_marker_test = TRUE
      i = outIndex
    }else{
      is_marker_test = FALSE
      stop("No markers are left in tractor_hybrid")
    }
  }else if(genoType == "bgen"){
"""
    text = replace_once_if_missing(
        text,
        insert_after_vcf,
        hybrid_branch,
        'genoType == "tractor_hybrid"',
        "SAIGE_SPATest_Marker.R hybrid setup",
    )
    text = text.replace('if(genoType == "plink"){\t\n      tempList = genoIndexList[[i]]', 'if(genoType == "plink"){\t\n      tempList = genoIndexList[[i]]', 1)
    text = text.replace(
        'if(genoType != "vcf"){\n      cat(paste0("(",Sys.time(),") ---- Analyzing Chunk ", i, "/", nChunks, ": chrom ", chrom," ---- \\\\n"))',
        'if(genoType != "vcf" && genoType != "tractor_hybrid"){\n      cat(paste0("(",Sys.time(),") ---- Analyzing Chunk ", i, "/", nChunks, ": chrom ", chrom," ---- \\\\n"))',
        1,
    )
    text = text.replace(
        'if(genoType != "vcf"){\n      cat(paste0("(",Sys.time(),") ---- Analyzing Chunk ',
        'if(genoType != "vcf" && genoType != "tractor_hybrid"){\n      cat(paste0("(",Sys.time(),") ---- Analyzing Chunk ',
        1,
    )
    text = must_replace(
        text,
        "  mainMarkerInCPP(genoType, traitType, genoIndex_prev, genoIndex, isMoreOutput, isImputation, isFirth)\n",
        """  if(genoType == "tractor_hybrid"){
    .tractor_io_before = get_TRACTORHYBRID_io_seconds()
    .tractor_decode_before = get_TRACTORHYBRID_decode_seconds()
    .tractor_get_marker_before = get_TRACTORHYBRID_get_marker_seconds()
    .tractor_pvalue_before = get_TRACTORHYBRID_pvalue_seconds()
    .tractor_output_before = get_TRACTORHYBRID_output_seconds()
    .tractor_reset_before = get_TRACTORHYBRID_reset_seconds()
    .tractor_condition_before = get_TRACTORHYBRID_condition_seconds()
    .tractor_impute_qc_before = get_TRACTORHYBRID_impute_qc_seconds()
    .tractor_joint_before = get_TRACTORHYBRID_joint_seconds()
    .tractor_variance_before = get_TRACTORHYBRID_variance_seconds()
    .tractor_condition_cache_hits_before = get_TRACTORHYBRID_condition_cache_hits()
    .tractor_condition_cache_misses_before = get_TRACTORHYBRID_condition_cache_misses()
    .tractor_saige_score_before = get_TRACTORHYBRID_saige_score_seconds()
    .tractor_saige_score_fast_calls_before = get_TRACTORHYBRID_saige_score_fast_calls()
    .tractor_saige_score_slow_calls_before = get_TRACTORHYBRID_saige_score_slow_calls()
    .tractor_saige_scorefast_extract_before = get_TRACTORHYBRID_saige_scorefast_extract_seconds()
    .tractor_saige_scorefast_projection_before = get_TRACTORHYBRID_saige_scorefast_projection_seconds()
    .tractor_saige_scorefast_variance_before = get_TRACTORHYBRID_saige_scorefast_variance_seconds()
    .tractor_saige_scorefast_result_before = get_TRACTORHYBRID_saige_scorefast_result_seconds()
    .tractor_saige_scorefast_gtilde_before = get_TRACTORHYBRID_saige_scorefast_gtilde_seconds()
    .tractor_saige_scorefast_gtilde_calls_before = get_TRACTORHYBRID_saige_scorefast_gtilde_calls()
    .tractor_saige_alloc_before = get_TRACTORHYBRID_saige_alloc_seconds()
    .tractor_saige_getadjg_before = get_TRACTORHYBRID_saige_getadjg_seconds()
    .tractor_saige_getadjg_calls_before = get_TRACTORHYBRID_saige_getadjg_calls()
    .tractor_saige_getadjg_accum_before = get_TRACTORHYBRID_saige_getadjg_accum_seconds()
    .tractor_saige_getadjg_projection_before = get_TRACTORHYBRID_saige_getadjg_projection_seconds()
    .tractor_saige_getadjg_spa_before = get_TRACTORHYBRID_saige_getadjg_spa_seconds()
    .tractor_saige_getadjg_spa_calls_before = get_TRACTORHYBRID_saige_getadjg_spa_calls()
    .tractor_saige_getadjg_firth_before = get_TRACTORHYBRID_saige_getadjg_firth_seconds()
    .tractor_saige_getadjg_firth_calls_before = get_TRACTORHYBRID_saige_getadjg_firth_calls()
    .tractor_saige_getadjg_condition_before = get_TRACTORHYBRID_saige_getadjg_condition_seconds()
    .tractor_saige_getadjg_condition_calls_before = get_TRACTORHYBRID_saige_getadjg_condition_calls()
    .tractor_saige_getadjg_region_before = get_TRACTORHYBRID_saige_getadjg_region_seconds()
    .tractor_saige_getadjg_region_calls_before = get_TRACTORHYBRID_saige_getadjg_region_calls()
    .tractor_saige_getadjg_other_before = get_TRACTORHYBRID_saige_getadjg_other_seconds()
    .tractor_saige_getadjg_other_calls_before = get_TRACTORHYBRID_saige_getadjg_other_calls()
    .tractor_saige_spa_before = get_TRACTORHYBRID_saige_spa_seconds()
    .tractor_saige_spa_calls_before = get_TRACTORHYBRID_saige_spa_calls()
    .tractor_saige_firth_before = get_TRACTORHYBRID_saige_firth_seconds()
    .tractor_saige_firth_calls_before = get_TRACTORHYBRID_saige_firth_calls()
    .tractor_saige_er_before = get_TRACTORHYBRID_saige_er_seconds()
    .tractor_saige_condition_before = get_TRACTORHYBRID_saige_condition_seconds()
    .tractor_saige_region_before = get_TRACTORHYBRID_saige_region_seconds()
    .tractor_chunk_time = system.time({
      mainMarkerInCPP(genoType, traitType, genoIndex_prev, genoIndex, isMoreOutput, isImputation, isFirth)
    })
    .tractor_io_delta = get_TRACTORHYBRID_io_seconds() - .tractor_io_before
    .tractor_decode_delta = get_TRACTORHYBRID_decode_seconds() - .tractor_decode_before
    .tractor_get_marker_delta = get_TRACTORHYBRID_get_marker_seconds() - .tractor_get_marker_before
    .tractor_pvalue_delta = get_TRACTORHYBRID_pvalue_seconds() - .tractor_pvalue_before
    .tractor_output_delta = get_TRACTORHYBRID_output_seconds() - .tractor_output_before
    .tractor_reset_delta = get_TRACTORHYBRID_reset_seconds() - .tractor_reset_before
    .tractor_condition_delta = get_TRACTORHYBRID_condition_seconds() - .tractor_condition_before
    .tractor_impute_qc_delta = get_TRACTORHYBRID_impute_qc_seconds() - .tractor_impute_qc_before
    .tractor_joint_delta = get_TRACTORHYBRID_joint_seconds() - .tractor_joint_before
    .tractor_variance_delta = get_TRACTORHYBRID_variance_seconds() - .tractor_variance_before
    .tractor_condition_cache_hits_delta = get_TRACTORHYBRID_condition_cache_hits() - .tractor_condition_cache_hits_before
    .tractor_condition_cache_misses_delta = get_TRACTORHYBRID_condition_cache_misses() - .tractor_condition_cache_misses_before
    .tractor_saige_score_delta = get_TRACTORHYBRID_saige_score_seconds() - .tractor_saige_score_before
    .tractor_saige_score_fast_calls_delta = get_TRACTORHYBRID_saige_score_fast_calls() - .tractor_saige_score_fast_calls_before
    .tractor_saige_score_slow_calls_delta = get_TRACTORHYBRID_saige_score_slow_calls() - .tractor_saige_score_slow_calls_before
    .tractor_saige_scorefast_extract_delta = get_TRACTORHYBRID_saige_scorefast_extract_seconds() - .tractor_saige_scorefast_extract_before
    .tractor_saige_scorefast_projection_delta = get_TRACTORHYBRID_saige_scorefast_projection_seconds() - .tractor_saige_scorefast_projection_before
    .tractor_saige_scorefast_variance_delta = get_TRACTORHYBRID_saige_scorefast_variance_seconds() - .tractor_saige_scorefast_variance_before
    .tractor_saige_scorefast_result_delta = get_TRACTORHYBRID_saige_scorefast_result_seconds() - .tractor_saige_scorefast_result_before
    .tractor_saige_scorefast_gtilde_delta = get_TRACTORHYBRID_saige_scorefast_gtilde_seconds() - .tractor_saige_scorefast_gtilde_before
    .tractor_saige_scorefast_gtilde_calls_delta = get_TRACTORHYBRID_saige_scorefast_gtilde_calls() - .tractor_saige_scorefast_gtilde_calls_before
    .tractor_saige_alloc_delta = get_TRACTORHYBRID_saige_alloc_seconds() - .tractor_saige_alloc_before
    .tractor_saige_getadjg_delta = get_TRACTORHYBRID_saige_getadjg_seconds() - .tractor_saige_getadjg_before
    .tractor_saige_getadjg_calls_delta = get_TRACTORHYBRID_saige_getadjg_calls() - .tractor_saige_getadjg_calls_before
    .tractor_saige_getadjg_accum_delta = get_TRACTORHYBRID_saige_getadjg_accum_seconds() - .tractor_saige_getadjg_accum_before
    .tractor_saige_getadjg_projection_delta = get_TRACTORHYBRID_saige_getadjg_projection_seconds() - .tractor_saige_getadjg_projection_before
    .tractor_saige_getadjg_spa_delta = get_TRACTORHYBRID_saige_getadjg_spa_seconds() - .tractor_saige_getadjg_spa_before
    .tractor_saige_getadjg_spa_calls_delta = get_TRACTORHYBRID_saige_getadjg_spa_calls() - .tractor_saige_getadjg_spa_calls_before
    .tractor_saige_getadjg_firth_delta = get_TRACTORHYBRID_saige_getadjg_firth_seconds() - .tractor_saige_getadjg_firth_before
    .tractor_saige_getadjg_firth_calls_delta = get_TRACTORHYBRID_saige_getadjg_firth_calls() - .tractor_saige_getadjg_firth_calls_before
    .tractor_saige_getadjg_condition_delta = get_TRACTORHYBRID_saige_getadjg_condition_seconds() - .tractor_saige_getadjg_condition_before
    .tractor_saige_getadjg_condition_calls_delta = get_TRACTORHYBRID_saige_getadjg_condition_calls() - .tractor_saige_getadjg_condition_calls_before
    .tractor_saige_getadjg_region_delta = get_TRACTORHYBRID_saige_getadjg_region_seconds() - .tractor_saige_getadjg_region_before
    .tractor_saige_getadjg_region_calls_delta = get_TRACTORHYBRID_saige_getadjg_region_calls() - .tractor_saige_getadjg_region_calls_before
    .tractor_saige_getadjg_other_delta = get_TRACTORHYBRID_saige_getadjg_other_seconds() - .tractor_saige_getadjg_other_before
    .tractor_saige_getadjg_other_calls_delta = get_TRACTORHYBRID_saige_getadjg_other_calls() - .tractor_saige_getadjg_other_calls_before
    .tractor_saige_spa_delta = get_TRACTORHYBRID_saige_spa_seconds() - .tractor_saige_spa_before
    .tractor_saige_spa_calls_delta = get_TRACTORHYBRID_saige_spa_calls() - .tractor_saige_spa_calls_before
    .tractor_saige_firth_delta = get_TRACTORHYBRID_saige_firth_seconds() - .tractor_saige_firth_before
    .tractor_saige_firth_calls_delta = get_TRACTORHYBRID_saige_firth_calls() - .tractor_saige_firth_calls_before
    .tractor_saige_er_delta = get_TRACTORHYBRID_saige_er_seconds() - .tractor_saige_er_before
    .tractor_saige_condition_delta = get_TRACTORHYBRID_saige_condition_seconds() - .tractor_saige_condition_before
    .tractor_saige_region_delta = get_TRACTORHYBRID_saige_region_seconds() - .tractor_saige_region_before
    .tractor_total = as.numeric(.tractor_chunk_time[["elapsed"]])
    .tractor_get_marker_overhead = max(0, .tractor_get_marker_delta - .tractor_io_delta - .tractor_decode_delta)
    .tractor_other = max(0, .tractor_total - .tractor_get_marker_delta - .tractor_pvalue_delta - .tractor_output_delta - .tractor_reset_delta - .tractor_impute_qc_delta - .tractor_joint_delta - .tractor_variance_delta)
    cat(sprintf("tractor_hybrid timing chunk %s: total=%.3fs input_io=%.3fs reader_decode=%.3fs hybrid_get_marker=%.3fs get_marker_overhead=%.3fs marker_pvalue=%.3fs reset_zero=%.3fs condition_total=%.3fs condition_cache_hit=%.0f condition_cache_miss=%.0f impute_qc=%.3fs variance_ratio=%.3fs joint_cct_spa=%.3fs output_write=%.3fs other_cpp_or_r=%.3fs\\n",
                i, .tractor_total, .tractor_io_delta, .tractor_decode_delta, .tractor_get_marker_delta, .tractor_get_marker_overhead, .tractor_pvalue_delta, .tractor_reset_delta, .tractor_condition_delta, .tractor_condition_cache_hits_delta, .tractor_condition_cache_misses_delta, .tractor_impute_qc_delta, .tractor_variance_delta, .tractor_joint_delta, .tractor_output_delta, .tractor_other))
    cat(sprintf("tractor_hybrid pvalue timing chunk %s: pvalue_score=%.3fs score_fast_calls=%.0f score_slow_calls=%.0f scorefast_extract=%.3fs scorefast_projection=%.3fs scorefast_variance=%.3fs scorefast_result=%.3fs scorefast_gtilde=%.3fs scorefast_gtilde_calls=%.0f pvalue_alloc=%.3fs pvalue_getadjg=%.3fs getadjg_calls=%.0f getadjg_accum=%.3fs getadjg_projection=%.3fs getadjg_spa=%.3fs getadjg_spa_calls=%.0f getadjg_firth=%.3fs getadjg_firth_calls=%.0f getadjg_condition=%.3fs getadjg_condition_calls=%.0f getadjg_region=%.3fs getadjg_region_calls=%.0f getadjg_other=%.3fs getadjg_other_calls=%.0f pvalue_spa=%.3fs spa_calls=%.0f pvalue_firth=%.3fs firth_calls=%.0f pvalue_er=%.3fs pvalue_condition_adjust=%.3fs pvalue_region_finalize=%.3fs\\n",
                i, .tractor_saige_score_delta, .tractor_saige_score_fast_calls_delta, .tractor_saige_score_slow_calls_delta, .tractor_saige_scorefast_extract_delta, .tractor_saige_scorefast_projection_delta, .tractor_saige_scorefast_variance_delta, .tractor_saige_scorefast_result_delta, .tractor_saige_scorefast_gtilde_delta, .tractor_saige_scorefast_gtilde_calls_delta, .tractor_saige_alloc_delta, .tractor_saige_getadjg_delta, .tractor_saige_getadjg_calls_delta, .tractor_saige_getadjg_accum_delta, .tractor_saige_getadjg_projection_delta, .tractor_saige_getadjg_spa_delta, .tractor_saige_getadjg_spa_calls_delta, .tractor_saige_getadjg_firth_delta, .tractor_saige_getadjg_firth_calls_delta, .tractor_saige_getadjg_condition_delta, .tractor_saige_getadjg_condition_calls_delta, .tractor_saige_getadjg_region_delta, .tractor_saige_getadjg_region_calls_delta, .tractor_saige_getadjg_other_delta, .tractor_saige_getadjg_other_calls_delta, .tractor_saige_spa_delta, .tractor_saige_spa_calls_delta, .tractor_saige_firth_delta, .tractor_saige_firth_calls_delta, .tractor_saige_er_delta, .tractor_saige_condition_delta, .tractor_saige_region_delta))
  }else{
    mainMarkerInCPP(genoType, traitType, genoIndex_prev, genoIndex, isMoreOutput, isImputation, isFirth)
  }
""",
        "SAIGE_SPATest_Marker.R timing wrapper",
        1,
    )
    text = must_replace(
        text,
        """  if(genoType == "vcf"){
    isEnd_Output =  check_Vcf_end()
  }else{
    isEnd_Output = (i==nChunks)\t
  }
""",
        """  if(genoType == "vcf"){
    isEnd_Output =  check_Vcf_end()
  }else if(genoType == "tractor_hybrid"){
    isEnd_Output = check_TRACTORHYBRID_end()
  }else{
    isEnd_Output = (i==nChunks)\t
  }
""",
        "SAIGE_SPATest_Marker.R end output",
        1,
    )
    text = must_replace(
        text,
        """  if(genoType == "vcf"){
    isVcfEnd =  check_Vcf_end()
    cat("isVcfEnd ", isVcfEnd, "\\n")
    if(isVcfEnd){
\tis_marker_test = FALSE\t     
    }
  }else{
    if(i > nChunks){
      is_marker_test = FALSE
    }\t    
  }
""",
        """  if(genoType == "vcf"){
    isVcfEnd =  check_Vcf_end()
    cat("isVcfEnd ", isVcfEnd, "\\n")
    if(isVcfEnd){
\tis_marker_test = FALSE\t     
    }
  }else if(genoType == "tractor_hybrid"){
    isHybridEnd = check_TRACTORHYBRID_end()
    cat("isTRACTORHYBRIDEnd ", isHybridEnd, "\\n")
    if(isHybridEnd){
      is_marker_test = FALSE
    }
  }else{
    if(i > nChunks){
      is_marker_test = FALSE
    }\t    
  }
""",
        "SAIGE_SPATest_Marker.R continuation",
        1,
    )
    write(path, text)


def patch_step2_cli():
    path = Path("/usr/local/bin/step2_SPAtests.R")
    text = read(path)
    text = replace_once_if_missing(
        text,
        """  make_option("--bgenFile", type="character",default="",
    help="Path to bgen file. Path to bgen file. Currently version 1.2 with 8 bit compression is supported"),
""",
        """  make_option("--bgenFile", type="character",default="",
    help="Path to bgen file. Path to bgen file. Currently version 1.2 with 8 bit compression is supported"),
  make_option("--tractorHybridPrefix", type="character",default="",
    help="Prefix for tractor_hybrid packed step2 input."),
""",
        "tractorHybridPrefix",
        "step2 cli option",
    )
    old = "bgenFile=opt$bgenFile,\n             bgenFileIndex=opt$bgenFileIndex,"
    new = old + "\n             tractorHybridPrefix=opt$tractorHybridPrefix,"

    first = text.find(old)
    if first < 0:
        raise RuntimeError("step2 cli main SPAGMMATtest call not found")
    if "tractorHybridPrefix=opt$tractorHybridPrefix" not in text[first:first + 300]:
        text = text[:first] + new + text[first + len(old):]

    second = text.find(old, first + len(new))
    if second >= 0 and "tractorHybridPrefix=opt$tractorHybridPrefix" not in text[second:second + 300]:
        text = text[:second] + new + text[second + len(old):]

    text = text.replace(
        "tractorHybridPrefix=opt$tractorHybridPrefix,\n             tractorHybridPrefix=opt$tractorHybridPrefix,",
        "tractorHybridPrefix=opt$tractorHybridPrefix,",
    )
    write(path, text)


def main():
    patch_main_hpp()
    patch_makevars()
    patch_main_cpp()
    patch_saige_test_cpp()
    patch_geno_r()
    patch_saige_test_main()
    patch_spatest_tractor()
    patch_spatest_marker()
    patch_step2_cli()


if __name__ == "__main__":
    main()
