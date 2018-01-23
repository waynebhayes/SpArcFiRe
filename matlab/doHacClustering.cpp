#include <limits.h>
#include <time.h>

#include <algorithm>
#include <list>
#include <string>
#include <vector>

#include "mex.h"
#include "matrix.h"

using std::list;
using std::sort;
using std::string;
using std::vector;

const int kMaxClusStartNbrs = 8;

// Describes a pixel-to-pixel similarity.
class SimlElt {
 public:
  SimlElt(int from, int to, double siml)
      : fromIdx(from), toIdx(to), simlVal(siml), valid(true) {
  }
        
  int getFromIdx() const {
    return fromIdx;
  }
  
  int getToIdx() const {
    return toIdx;
  }
  
  double getSimlVal() const {
    return simlVal;
  }
  
  // Determines whether this pixel similarity pair has been invalidated by
  // a previous merge rejection; i.e. (1) it would have been removed by 
  // similarity-matrix-based HAC, (2) the corresponding clusters' merge was 
  // rejected by the arc-fit criterion, and (3) the similarity linked the
  // two clusters at the time of the rejection. A stricter definition of 
  // "invalid" would include only part (1), but most of the time we can 
  // simply check to see if the similarity is now within the same cluster
  // before merging. By explicitly tracking invalid merges under all three
  // criteria, we can avoid arc-merge checks that would not have happened
  // under similarity matrix-based HAC, maintaining result-equivalence and 
  // saving time.
  bool IsValid() const {
    return valid;
  }
  
  // Marks this similarity as invalid; see IsInvalid() for more information.
  void SetInvalid() {
    valid = false;
  }
  
  bool equals(const SimlElt &other) const {
    return (this->fromIdx == other.fromIdx) && 
      (this->toIdx == other.toIdx) && (this->simlVal == other.simlVal) &&
      (this->valid == other.valid);
  }
  
  void display() const {
//    char *buf[1024];
    mexPrintf("SimlElt fromIdx=%d toIdx=%d simlVal=%.4f\n", fromIdx, 
        toIdx, simlVal);
  }
  
 private:
  int fromIdx;
  int toIdx;
  double simlVal;
  bool valid;
};

class SimlEltComparator {
 public:
  // Determines which similarity should be viewed first (i.e., which 
  // similarity is larger). Note that this is the opposite of what would 
  // be expected from operator< .
  bool operator()(const SimlElt &se1, const SimlElt &se2) const {
//    if ((se1 == NULL) || (se2 == NULL)) {
//      mexPrintf("internal error: null SimlElt\n");
//    }
    return (se1.getSimlVal() > se2.getSimlVal());
  }
};

class SimlPosnComparator {
 public:
  // main_siml_list must not be destroyed before SimlPosnComparator
  SimlPosnComparator(
      const vector<SimlElt> &main_siml_list) 
      : main_siml_list_(main_siml_list) {
  }
  // Determines which similarity should be viewed first (i.e., which 
  // similarity is larger). Note that this is the opposite of what would 
  // be expected from operator< .
  bool operator()(const size_t idx1, 
                                      const size_t idx2) const {
    return (main_siml_list_.at(idx1).getSimlVal() > 
        main_siml_list_.at(idx2).getSimlVal());
  }
 private:
  const vector<SimlElt> &main_siml_list_;
};

class ImgPxl {
 public:
  ImgPxl(mwIndex img_idx, double brt) 
      : img_idx_(img_idx), brt_(brt) {
    mxAssert(brt != 0, "cluster pixel brightness should not be zero");
  }
  
  mwIndex getIdx() const {
    return img_idx_;
  }
          
  double getBrt() const {
    return brt_;
  }

 private:
  mwIndex img_idx_;
  double brt_;
};

// Tracks cluster pixels
// TODO: extend to calculate and remember log-spiral fits
class HacCluster {
 public:
  HacCluster(ImgPxl startPxl)
      : params_(NULL), bounds_(NULL) {
    pxls_ = vector<ImgPxl>(1, startPxl);
  }
  
  // Adds the other cluster's contents to this cluster. Does not free the 
  // memory for the other cluster; the caller is responsible for this.
  void absorb(const HacCluster &other) {
    pxls_.insert(pxls_.end(), other.pxls_.begin(), other.pxls_.end());
    // TODO: delete params and bounds once functionality 
    // is added to manage these
  }
  
  // Create a cluster representation as a mxArray of doubles. The caller 
  // is responsible for freeing the allocated memory.
  mxArray *genClusMtx(mwSize num_rows, mwSize num_columns) 
      const {
    mxArray *clus_mtx = mxCreateDoubleMatrix(num_rows, num_columns, 
                                             mxREAL);
    double *pr = static_cast<double*>(mxGetData(clus_mtx));
    mxAssert(!pxls_.empty(), "pixel cluster should not be empty");
    for (vector<ImgPxl>::const_iterator it = pxls_.begin(); 
        it < pxls_.end(); ++it) {
      ImgPxl curPxl = *it;
      *(pr + curPxl.getIdx() - 1) = curPxl.getBrt();
      mxAssert(curPxl.getBrt() != 0, "pixel brightness should not be zero");
    }
    return clus_mtx;
  }
  
  mwSize numPts() const {
    return pxls_.size();
  }

 private:
  vector<ImgPxl> pxls_;
  mxArray *params_;
  mxArray *bounds_;
};

void chkArgsNonNull(int nrhs, const mxArray *prhs[]) {
  for (mwSize i=0; i<nrhs; ++i) {
    if (prhs[i] == NULL) {
      string msg("mexFunction inputs cannot be null");
      mexErrMsgTxt(msg.c_str());
    }
  }
}

void chkIsRealFullNumeric(const mxArray *var, const string *varname) {
  if (mxIsComplex(var)) {
    string msg("input '" + *varname + "' must be real-valued");
    mexErrMsgTxt(msg.c_str());
  }
  if (mxIsSparse(var)) {
    string msg("input '" + *varname + "' must be non-sparse");
    mexErrMsgTxt(msg.c_str());
  }
  if (!mxIsNumeric(var)) {
    string msg("input '" + *varname + "' must be numeric");
    mexErrMsgTxt(msg.c_str());
  }
}

void chkIsStruct(const mxArray *var, const string *varname) {
  if (!mxIsStruct(var)) {
    string msg("input '" + *varname + "' must be a struct");
    mexErrMsgTxt(msg.c_str());
  }
}

void chkSameDims(const mxArray *arr1, const mxArray *arr2, 
    const string *varname1, const string *varname2) {
  mwSize nDims = mxGetNumberOfDimensions(arr1);
  if (nDims != mxGetNumberOfDimensions(arr2)) {
    string msg("inputs '" + *varname1 + "' and '" + *varname2 + 
        "' have a different number of dimensions");
    mexErrMsgTxt(msg.c_str());
  }
  
  const mwSize *dims1 = mxGetDimensions(arr1);
  const mwSize *dims2 = mxGetDimensions(arr2);
  for (mwSize i=0; i<nDims; ++i) {
    if (dims1[i] != dims2[i]) {
      string msg("inputs '" + *varname1 + "' and '" + *varname2 + 
          "' have a different size");
      mexErrMsgTxt(msg.c_str());
    }
  }
}

// Manages assignments of pixels to clusters, calculates arc-fit cost of
// merging two pairs of clusters, and tracks invalidated entries on the 
// merge list
class ClusPartition {
 public:
  // Creates the cluster partition by assigning a singleton cluster to
  // each pixel (that has at least one nonzero similarity, as indicated by 
  // add_as_cluster).
  // main_siml_list must not be destroyed before ClusPartition
  ClusPartition(const mxArray *img, 
                const vector<bool> &add_as_cluster, 
                vector<SimlElt> &main_siml_list,
                const mxArray *ctr_row,
                const mxArray *ctr_column,
                const mxArray *bar_info,
                const mxArray *stgs) 
      : asgns_(NULL), main_siml_list_(main_siml_list) {
    Init(img, add_as_cluster, main_siml_list, ctr_row, ctr_column, 
        bar_info, stgs);
  }
 
  ~ClusPartition() {
    mxDestroyArray(ctr_row_);
    mxDestroyArray(ctr_column_);
    mxDestroyArray(bar_info_);
    mxDestroyArray(stgs_);
    
    for (vector<HacCluster*>::iterator iter = clusts_.begin();
        iter < clusts_.end(); ++iter) {
      if (*iter != NULL) {
        delete *iter;
      }
    }
    
    for (vector<list<size_t>*>::iterator iter = perClusterSimlPosns_.begin();
        iter < perClusterSimlPosns_.end(); ++iter) {
      if (*iter != NULL) {
        delete *iter;
      }
    }
    
    delete simlCmp_;
  }
  
  // Determines whether the two given cluster representatives belong to 
  // the same cluster
  bool HasSameCluster(mwIndex rep1, mwIndex rep2) {
    mwIndex asgn1 = GetAssignment(rep1);
    mwIndex asgn2 = GetAssignment(rep2);
    return (asgn1 == asgn2);
  }
  
  // Merges two clusters containing pixels rep1 and rep2. These can be 
  // indexes of any pixel in the cluster.
  //
  // Returns true if a merge was performed (i.e., the two members were 
  // from unique clusters.
  bool Merge(const SimlElt &mergeSiml) {
    mwIndex rep1 = mergeSiml.getFromIdx();
    mwIndex rep2 = mergeSiml.getToIdx();
    if (HasSameCluster(rep1, rep2)) {
      SkipMerge(mergeSiml);
      return false;
    }
    mwIndex asgn1 = GetAssignment(rep1);
    mwIndex asgn2 = GetAssignment(rep2);
    HacCluster *clus1 = clusts_.at(asgn1);
    HacCluster *clus2 = clusts_.at(asgn2);
    if ((clus1 == NULL) || (clus2 == NULL)) {
      string msg("internal error: null cluster object at cluster ");
      msg.append("assignment index");
      mexErrMsgTxt(msg.c_str());
    }
    
    list<size_t> *posns1 = perClusterSimlPosns_.at(asgn1);
    list<size_t> *posns2 = perClusterSimlPosns_.at(asgn2);
    if ((posns1 == NULL) || (posns2 == NULL)) {
      char msg[256];
      sprintf(msg, "Internal Error: null posns1 or posns2 (%p,\t%p)",
          posns1, posns2);
      mexErrMsgTxt(msg);
    }
    if (!main_siml_list_.at(posns1->front()).equals(mergeSiml)) {
      main_siml_list_.at(posns1->front()).display();
      main_siml_list_.at(posns2->front()).display();
      mergeSiml.display();
      mexErrMsgTxt(
          "Internal Error: main-list next siml differs from per-cluster next siml");
    } 
    posns1->pop_front();
    if (clus1->numPts() > clus2->numPts()) {
      clus1->absorb(*clus2);
      delete clus2;
      clusts_.at(asgn2) = NULL;
      asgns_.at(asgn2) = asgn1;
      posns1->merge(*posns2, *simlCmp_);
      delete perClusterSimlPosns_.at(asgn2);
      perClusterSimlPosns_.at(asgn2) = NULL;
    } else {
      clus2->absorb(*clus1);
      delete clus1;
      clusts_.at(asgn1) = NULL;
      asgns_.at(asgn1) = asgn2;
      posns2->merge(*posns1, *simlCmp_);
      delete perClusterSimlPosns_.at(asgn1);
      perClusterSimlPosns_.at(asgn1) = NULL;
    }
  }
  
  // Updates the cluster's similarity list to reflect that the next 
  // similarity was skipped (since the points are now in the same cluster, 
  // or the merge was rejected)
  void SkipMerge(const SimlElt &mergeSiml) {
    mwIndex rep1 = mergeSiml.getFromIdx();
    mwIndex rep2 = mergeSiml.getToIdx();
    mwIndex asgn1 = GetAssignment(rep1);
    mwIndex asgn2 = GetAssignment(rep2);
    
    list<size_t> *posns1 = perClusterSimlPosns_.at(asgn1);
    list<size_t> *posns2 = perClusterSimlPosns_.at(asgn2);
    if ((posns1 == NULL) || (posns2 == NULL)) {
      char msg[256];
      sprintf(msg, "Internal Error: null posns1 or posns2 (%p,\t%p)",
          posns1, posns2);
      mexErrMsgTxt(msg);
    }
    if (!mergeSiml.IsValid()) {
      // invalidated similarity already deleted from cluster list
      // no need to do anything
      // Note: this case must be handled first, since posns1 could be empty
    } else if (main_siml_list_.at(posns1->front()).equals(mergeSiml)) {
      // keep the per-cluster similarities updated by removing the
      // skipped similarity
      posns1->pop_front();
    } else {
      main_siml_list_.at(posns1->front()).display();
      main_siml_list_.at(posns2->front()).display();
      mergeSiml.display();
      mexErrMsgTxt(
          "Internal Error: main-list next siml differs from per-cluster next siml");
    } 
  }
  
  // Performs updates to account for a merge rejection between clusters
  // with representatives clus1rep and clus2rep (due to the arc-fit 
  // criterion). Similarities from one cluster to the other will be 
  // invalidated (they would be removed by the max operation in the 
  // straightforward HAC implementation and can be skipped when they point 
  // within the same cluster, but we want to do an explicit and earlier 
  // invalidation to avoid additional arc-fit merge checks.
  void RecordMergeRejection(mwIndex clus1rep, 
                                           mwIndex clus2rep) {
    RecordMergeRejectionOneDirectional(clus1rep, clus2rep);
    RecordMergeRejectionOneDirectional(clus2rep, clus1rep);
  }
  
  // Determines whether the given (pixel) index is assigned to a cluster.
  bool HasAssignment(mwIndex index) const {
    bool hasAsgn = ((index >= 1) && (index < asgns_.size()) && 
        (asgns_.at(index) >= 1));
    return hasAsgn;
  }
  
  // Determines the current cluster representative for the given 
  // (pixel) index
  mwIndex GetAssignment(mwIndex index) {
    // check that ind is in range; 0th index should never be used since 
    // we're aligning with Matlab's indexing
    if ((index < 1) || (index >= asgns_.size())) {
      string msg("internal error: out-of-range index ");
      char buf[50];
      sprintf(buf, "(%d)", index);
      msg.append(buf);
      msg.append(" when getting cluster assignment [asgns size = ");
      sprintf(buf, "%d]", asgns_.size());
      msg.append(buf);
      mexErrMsgTxt(msg.c_str());
    }
    if (asgns_.at(index) < 1) {
      string msg("internal error: attempted to get assignment for ");
      msg.append("unassigned pixel index ");
      char buf[50];
      sprintf(buf, "(%d)", index);
      msg.append(buf);
      mexErrMsgTxt(msg.c_str());
    }
    
    mwIndex prev_index = index;
    mwIndex cur_index = asgns_.at(prev_index);
    mwIndex next_index = asgns_.at(cur_index);
    // avoid vector allocation if we don't need to trace our path; looking 
    // one index ahead so we can do this for all members that point 
    // directly to the final index, not just the pixel that also has the 
    // final index
    if (cur_index == next_index) {
      // then asgns_(index) gives final assignment
      return cur_index;
    }
    vector<mwIndex> indirect_indexes = vector<mwIndex>(1, prev_index);
    while (cur_index != next_index) {
      indirect_indexes.push_back(cur_index);
      cur_index = next_index;
      next_index = asgns_.at(cur_index);
    }
    if (cur_index < 1) {
      mexErrMsgTxt("internal error: update to unassigned index");
    }
    // union-find path compression
    vector<mwIndex>::iterator it;
    for (it = indirect_indexes.begin(); it < indirect_indexes.end(); 
        ++it) {
      asgns_.at(*it) = cur_index;
    }
    return cur_index;
  }
  
  mwSize GetClusterSize(mwIndex clus_rep_index) {
    mwIndex asgn_idx = GetAssignment(clus_rep_index);
    return (clusts_.at(asgn_idx))->numPts();
  }
  
  // Determines the arc-fit merge error
  double CalcArcMergeCost(mwIndex rep_index_1, 
      mwIndex rep_index_2) {
    mwIndex asgn_idx_1 = GetAssignment(rep_index_1);
    mwIndex asgn_idx_2 = GetAssignment(rep_index_2);
    mxArray *clus_mtx_1 = (clusts_.at(asgn_idx_1))->genClusMtx(
        num_image_rows_, num_image_columns_);
    mxArray *clus_mtx_2 = (clusts_.at(asgn_idx_2))->genClusMtx(
        num_image_rows_, num_image_columns_);
    
    mxArray *lhs[1];
    mxArray *rhs[12];
    mxArray *empty = mxCreateDoubleMatrix(0, 0, mxREAL);
    rhs[0] = clus_mtx_1;
    rhs[1] = rhs[2] = rhs[3] = empty;
    rhs[4] = clus_mtx_2;
    rhs[5] = rhs[6] = rhs[7] = empty;
    rhs[8] = ctr_row_;
    rhs[9] = ctr_column_;
    rhs[10] = bar_info_;
    rhs[11] = stgs_;
//    time_t start,end;
//    time(&start);
    int call_status = mexCallMATLAB(1, lhs, 12, rhs, "calcArcMergeErr");
//    time(&end);
//    mexPrintf("time for mexCallMatlab for merge error calculation: %f\n",
//        difftime(end, start));
    mxDestroyArray(empty);
    if (call_status != 0) {
      char buf[256];
      sprintf(buf, "call to 'calcArcMergeErr' failed (code %d)", 
          call_status);
      mexErrMsgTxt(buf);
    }
    
    double *p_cost = static_cast<double*>(mxGetData(lhs[0]));
    
    mxDestroyArray(clus_mtx_1);
    mxDestroyArray(clus_mtx_2);
    
    return *p_cost;
  }
  
  // create a mxArray, with the same dimensions as the input image,
  // where each element corresponds to the pixel's cluster assignment 
  // (zero if no such assignment). The caller is responsible for destroying 
  // the mxArray.
  mxArray *buildClusterAssignmentMatrix() {
    mxArray *assignment_matrix = mxCreateNumericMatrix(num_image_rows_, 
        num_image_columns_, mxINT32_CLASS, mxREAL);
    int *pr = static_cast<int*>(mxGetData(assignment_matrix));
    
    mwSize num_img_elements = num_image_rows_ * num_image_columns_;
    for (mwIndex i=1; i <= num_img_elements; ++i) {
      if (HasAssignment(i)) {
        *pr = GetAssignment(i);
      }
      ++pr;
    }
    return assignment_matrix;
  }
  
 private:
  void Init(const mxArray *img,
                           const vector<bool> &add_as_cluster,
                           const vector<SimlElt> &main_siml_list,
                           const mxArray *ctr_row,
                           const mxArray *ctr_column,
                           const mxArray *bar_info,
                           const mxArray *stgs) {
    if (!mxIsDouble(img)) {
      mexErrMsgTxt("'img' values must be of type double");
    }
    if (mxGetNumberOfDimensions(img) != 2) {
      mexErrMsgTxt("'img' must be 2D array");
    }
    
    num_image_rows_ = mxGetM(img);
    num_image_columns_ = mxGetN(img);
    if ((num_image_rows_ * num_image_columns_) > LONG_MAX) {
      mexErrMsgTxt("image size greater than LONG_MAX on this platform\n");
    }
//    ctr_row_ = mxCreateDoubleMatrix(1, 1, mxREAL);
//    double *pr = static_cast<double*>(mxGetData(ctr_row_));
//    *pr = num_image_rows_/2;
//    ctr_column_ = mxCreateDoubleMatrix(1, 1, mxREAL);
//    pr = static_cast<double*>(mxGetData(ctr_column_));
//    *pr = num_image_columns_/2;
    ctr_row_ = mxDuplicateArray(ctr_row);
    ctr_column_ = mxDuplicateArray(ctr_column);
    
    mwSize nElt = mxGetNumberOfElements(img);
    asgns_ = vector<long>(nElt+1, -1);
    clusts_ = vector<HacCluster*>(nElt+1, NULL);
    
    // 1-indexed to correspond with Matlab
    perClusterSimlPosns_ = vector<list<size_t>*>(nElt+1, NULL);  
    for (size_t i = 0; i < main_siml_list.size(); ++i) {
      SimlElt curSiml = main_siml_list.at(i);
      int curFromIdx = curSiml.getFromIdx();
      if (perClusterSimlPosns_.at(curFromIdx) == NULL) {
        perClusterSimlPosns_.at(curFromIdx) = 
            new list<size_t>();
      }
      // since main_siml_list is sorted, individual cluster lists 
      // will also be sorted
      perClusterSimlPosns_.at(curFromIdx)->push_back(i);
      // only one cluster needs to hold the similarity; we always make it
      // the cluster given in fromIdx
    }
    
    double *pr = static_cast<double*>(mxGetData(img));
//    pr = static_cast<double*>(mxGetData(img));
    for (mwIndex i=0; i<nElt; ++i) {
      double curImgVal = *pr;
      if (add_as_cluster.at(i+1)) {
        // Matlab uses 1-indexing; keep correspondence with indices 
        // given by similarity pairs
        asgns_.at(i+1) = i+1;
        clusts_.at(i+1) = 
            new HacCluster(ImgPxl(i+1, curImgVal));
        if (perClusterSimlPosns_.at(i+1) == NULL) {
          // cluster has nonzero points, but all its similarities are 
          // tracked by other clusters
          perClusterSimlPosns_.at(i+1) = new list<size_t>();
        }
      }
      ++pr;
    }
    
    bar_info_ = mxDuplicateArray(bar_info);
    stgs_ = mxDuplicateArray(stgs);
    simlCmp_ = new SimlPosnComparator(main_siml_list);
  }
  
  // see RecordMergeRejection()
  void RecordMergeRejectionOneDirectional(
      mwIndex clus1rep, 
      mwIndex clus2rep) {
    mwIndex asgn1 = GetAssignment(clus1rep);
    list<size_t> *simlPosns1 = perClusterSimlPosns_.at(asgn1);
    for (list<size_t>::iterator iter = simlPosns1->begin(); 
        iter != simlPosns1->end(); ++iter) {
      size_t cur_posn = *iter;
      SimlElt &cur_siml = main_siml_list_.at(cur_posn);
      int from_idx = cur_siml.getFromIdx();
      int to_idx = cur_siml.getToIdx();
      if (HasSameCluster(from_idx, to_idx)) {
        // avoid duplicate examination of within-cluster similarities 
        // during bad-merge rejections; can help running-time analysis
        cur_siml.SetInvalid();
        iter = simlPosns1->erase(iter);
        --iter;
      } else if (HasSameCluster(to_idx, clus2rep)) {
        // this similarity has been invalidated by the merge rejection
        cur_siml.SetInvalid();
        iter = simlPosns1->erase(iter);
        --iter;
      }
    }
  }
  
  mwSize num_image_rows_;
  mwSize num_image_columns_;
  mxArray *ctr_row_;
  mxArray *ctr_column_;
  vector<long> asgns_; 
  vector<HacCluster*> clusts_;
  mxArray *bar_info_;
  mxArray *stgs_;
  vector<list<size_t>*> perClusterSimlPosns_;
  vector<SimlElt> &main_siml_list_;
  
  SimlPosnComparator *simlCmp_;
};

// Assumes that rIdxs, cIdxs, and allSimls are all MATLAB vectors with 
// the same size
mxArray *doHacClustering(const mxArray *rIdxs, 
                         const mxArray *cIdxs, 
                         const mxArray *allSimls,
                         const mxArray *img,
                         const mxArray *ctrR,
                         const mxArray *ctrC,
                         const mxArray *barInfo,
                         const mxArray *stgs) { 
  // get setting "mergeChkMinClusSz" from Matlab stgs struct
  mxArray *arr_stg_val = mxGetField(stgs, 0, "mergeChkMinClusSz");
  if (arr_stg_val == NULL) {
    mexErrMsgTxt("could not find 'mergeChkMinClusSz' in 'stgs'");
  }
  double *p_stg_val = static_cast<double*>(mxGetData(arr_stg_val));
  int min_clus_size_for_merge_check = static_cast<int>(*p_stg_val);
  char buf[128];
  sprintf(buf, "min_clus_size_for_merge_check = %d\n", 
      min_clus_size_for_merge_check);
  mexPrintf(buf);
  
  // get setting "errRatioThres" from Matlab stgs struct
  arr_stg_val = mxGetField(stgs, 0, "errRatioThres");
  if (arr_stg_val == NULL) {
    mexErrMsgTxt("could not find 'errRatioThres' in 'stgs'");
  }
  p_stg_val = static_cast<double*>(mxGetData(arr_stg_val));
  double err_ratio_thres = *p_stg_val;
  char buf2[128];
  sprintf(buf2, "err_ratio_thres = %2.4f\n", err_ratio_thres);
  mexPrintf(buf2);
  
  mwSize nInSimls = mxGetNumberOfElements(rIdxs);
  vector<SimlElt> mainSimlList = vector<SimlElt>();
  mainSimlList.reserve(static_cast<size_t>(nInSimls));
 
  double *pRIdx = mxGetPr(rIdxs);
  double *pCIdx = mxGetPr(cIdxs);
  double *pSml = mxGetPr(allSimls);
  mwSize nPix = mxGetNumberOfElements(img);
  // Which pixels have at least one nonzero similarity associated with 
  // them, and should thus have a corresponding initial cluster.
  // 1-indexed to correspond with Matlab
  vector<bool> add_as_cluster(nPix+1, false);
  for (mwSize i=0; i < nInSimls; ++i) {
    int curRIdx = static_cast<int>(*pRIdx);
    int curCIdx = static_cast<int>(*pCIdx);
    double curSml = *pSml;
    if ((curRIdx < curCIdx) && (curSml > 0)) {
        // similarity is non-redundant and non-self, so add it
        mainSimlList.push_back(SimlElt(curRIdx, curCIdx, curSml));
        // pixels have at least one nonzero similarity associated with them
        add_as_cluster.at(curRIdx) = true;
        add_as_cluster.at(curCIdx) = true;
    }
    ++pRIdx;
    ++pCIdx;
    ++pSml;
  }
  
  sort(mainSimlList.begin(), mainSimlList.end(), SimlEltComparator());
  
  ClusPartition *partition = new ClusPartition(img, add_as_cluster, 
      mainSimlList, ctrR, ctrC, barInfo, stgs);
  
  mexPrintf("performing clustering...\n");
  for (vector<SimlElt>::iterator it=mainSimlList.begin(); 
      it < mainSimlList.end(); ++it) {
    SimlElt curSiml = *it;
    mwIndex fromIdx = curSiml.getFromIdx();
    mwIndex toIdx = curSiml.getToIdx();
    if (!curSiml.IsValid() || partition->HasSameCluster(fromIdx, toIdx)) {
      partition->SkipMerge(curSiml);
      continue;
    }
    mwSize fromSize = partition->GetClusterSize(fromIdx);
    mwSize toSize = partition->GetClusterSize(toIdx);
    bool above_min_size = (fromSize > min_clus_size_for_merge_check) &&
        (toSize > min_clus_size_for_merge_check);
    bool merge_cost_ok = true;
    if (above_min_size) {
      double merge_cost = 
          partition->CalcArcMergeCost(fromIdx, toIdx);
      merge_cost_ok = (merge_cost <= err_ratio_thres);
    }
    if (!above_min_size || merge_cost_ok) {
      partition->Merge(curSiml);
    } else {
      partition->SkipMerge(curSiml);
    }
    if (above_min_size && !merge_cost_ok) {
      partition->RecordMergeRejection(fromIdx, toIdx);
    } 
  }
  
  mxArray *assignment_matrix = partition->buildClusterAssignmentMatrix();
  delete partition;
  return assignment_matrix;
}

// input arguments should be: rIdxs, cIdxs, allSimls, img, ctrR, ctrC, barInfo, stgs
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  if (nrhs != 8) {
    char buf[50]; sprintf(buf, "%d", nrhs);
    string msg("doHacClustering.cpp: expected 6 arguments, received ");
    msg.append(buf);
    mexErrMsgTxt(msg.c_str());
  }
  
  if (nlhs != 1) {
    mexErrMsgTxt("doHacClustering.cpp: should have 1 output argument");
  }
  
  chkArgsNonNull(nrhs, prhs);
  
  string varname("rIdxs"); chkIsRealFullNumeric(prhs[0], &varname);
  varname = "cIdxs"; chkIsRealFullNumeric(prhs[1], &varname);
  varname = "allSimls"; chkIsRealFullNumeric(prhs[2], &varname);
  varname = "img"; chkIsRealFullNumeric(prhs[3], &varname);
  varname = "ctrR"; chkIsRealFullNumeric(prhs[4], &varname);
  varname = "ctrC"; chkIsRealFullNumeric(prhs[5], &varname);
  varname = "barInfo"; chkIsStruct(prhs[6], &varname);
  varname = "stgs"; chkIsStruct(prhs[7], &varname);
  
  string varname0("rIdxs"); 
  string varname1("cIdxs"); 
  chkSameDims(prhs[0], prhs[1], &varname0, &varname1);
  string varname2("allSimls");
  chkSameDims(prhs[1], prhs[2], &varname1, &varname2);
  
  plhs[0] = doHacClustering(prhs[0], prhs[1], prhs[2], prhs[3], prhs[4], 
      prhs[5], prhs[6], prhs[7]);
}