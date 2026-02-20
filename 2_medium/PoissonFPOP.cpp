// PoissonFPOP.cpp
// Unconstrained optimal partitioning (Poisson loss) via FPOP.
// Based on PeakSegOptimal's constrained solver -- simplified from
// 2-state (up/down) to 1-state by replacing the constrained min
// operators with a single unconstrained min.
// Piecewise function code from: github.com/tdhock/PeakSegOptimal

#include <Rcpp.h>
#include <list>
#include <math.h>
#include <vector>

#define NEWTON_EPSILON 1e-12
#define NEWTON_STEPS 100
#define PREV_NOT_SET (-3)
#define ABS(x) ((x)<0 ? -(x) : (x))
#define ERROR_MIN_MAX_SAME 1
#define ERROR_NEGATIVE_DATA 2

// --- data structures (from funPieceListLog.h) ---

// single piece: Linear*exp(x) + Log*x + Constant in log-mean space
class PoissonLossPieceLog {
public:
  double Linear;
  double Log;
  double Constant;
  double min_log_mean;
  double max_log_mean;
  int data_i;
  double prev_log_mean;
  PoissonLossPieceLog();
  PoissonLossPieceLog(double li, double lo, double co,
                      double m, double M, int i, double prev);
  double argmin();
  double argmin_mean();
  double getCost(double log_mean);
  double getDeriv(double log_mean);
  double PoissonLoss(double mean);
  double PoissonDeriv(double mean);
  bool has_two_roots(double equals);
  double get_smaller_root(double equals);
  double get_larger_root(double equals);
  void print();
};

typedef std::list<PoissonLossPieceLog> PoissonLossPieceListLog;

class PiecewisePoissonLossLog {
public:
  PoissonLossPieceListLog piece_list;
  void set_to_unconstrained_min_of(PiecewisePoissonLossLog *input);
  void set_to_min_env_of(PiecewisePoissonLossLog *fun1,
                         PiecewisePoissonLossLog *fun2, int verbose);
  void push_min_pieces(PiecewisePoissonLossLog *fun1,
                       PiecewisePoissonLossLog *fun2,
                       PoissonLossPieceListLog::iterator it1,
                       PoissonLossPieceListLog::iterator it2, int verbose);
  void push_piece(PoissonLossPieceListLog::iterator it,
                  double min_log_mean, double max_log_mean);
  void add(double Linear, double Log, double Constant);
  void multiply(double x);
  void set_prev_seg_end(int prev_seg_end);
  void findMean(double log_mean, int *seg_end, double *prev_log_mean);
  double findCost(double log_mean);
  void Minimize(double *best_cost, double *best_log_mean,
                int *data_i, double *prev_log_mean);
  void print();
};

bool sameFuns(PoissonLossPieceListLog::iterator it1,
              PoissonLossPieceListLog::iterator it2);

// --- PoissonLossPieceLog methods (unchanged from reference) ---

PoissonLossPieceLog::PoissonLossPieceLog
(double li, double lo, double co, double m, double M, int i, double prev){
  Linear = li; Log = lo; Constant = co;
  min_log_mean = m; max_log_mean = M;
  data_i = i; prev_log_mean = prev;
}

PoissonLossPieceLog::PoissonLossPieceLog(){}

double PoissonLossPieceLog::argmin_mean(){
  return -Log / Linear;
}

double PoissonLossPieceLog::argmin(){
  return log(argmin_mean());
}

double PoissonLossPieceLog::getCost(double log_mean){
  double linear_term, log_term;
  if(log_mean == INFINITY){
    return (0 < Linear) ? INFINITY : -INFINITY;
  }
  if(log_mean == -INFINITY){
    linear_term = 0.0;
  }else{
    linear_term = Linear*exp(log_mean);
  }
  if(Log==0){
    log_term = 0.0;
  }else{
    log_term = Log*log_mean;
  }
  return linear_term + log_term + Constant;
}

double PoissonLossPieceLog::getDeriv(double log_mean){
  double linear_term;
  if(log_mean == -INFINITY){
    linear_term = 0.0;
  }else{
    linear_term = Linear*exp(log_mean);
  }
  return linear_term + Log;
}

double PoissonLossPieceLog::PoissonLoss(double mean){
  double loss_without_log_term = Linear*mean + Constant;
  if(Log==0) return loss_without_log_term;
  return loss_without_log_term + Log*log(mean);
}

double PoissonLossPieceLog::PoissonDeriv(double mean){
  return Linear + Log/mean;
}

bool PoissonLossPieceLog::has_two_roots(double equals){
  if(Log == 0){
    throw "degenerate: Log==0 in has_two_roots";
    return false;
  }
  double optimal_mean = argmin_mean();
  double optimal_log_mean = log(optimal_mean);
  double optimal_cost = getCost(optimal_log_mean);
  double optimal_cost2 = PoissonLoss(optimal_mean);
  if(0 < Linear){
    return optimal_cost + NEWTON_EPSILON < equals &&
           optimal_cost2 + NEWTON_EPSILON < equals;
  }
  return equals + NEWTON_EPSILON < optimal_cost &&
         equals + NEWTON_EPSILON < optimal_cost2;
}

double PoissonLossPieceLog::get_larger_root(double equals){
  double optimal_mean = argmin_mean();
  double optimal_cost = PoissonLoss(optimal_mean);
  double right_cost = getCost(max_log_mean);
  if((optimal_cost < right_cost && right_cost < equals) ||
     (optimal_cost > right_cost && right_cost > equals)){
    return max_log_mean+1;
  }
  double candidate_root = optimal_mean + 1;
  double candidate_cost, possibly_outside, deriv;
  double closest_positive_cost = INFINITY, closest_positive_mean;
  double closest_negative_cost = -INFINITY, closest_negative_mean;
  if(optimal_cost < 0){
    closest_negative_cost = optimal_cost;
    closest_negative_mean = optimal_mean;
  }else{
    closest_positive_cost = optimal_cost;
    closest_positive_mean = optimal_mean;
  }
  int step=0;
  do{
    candidate_cost = PoissonLoss(candidate_root) - equals;
    if(0 < candidate_cost && candidate_cost < closest_positive_cost){
      closest_positive_cost = candidate_cost;
      closest_positive_mean = candidate_root;
    }
    if(closest_negative_cost < candidate_cost && candidate_cost < 0){
      closest_negative_cost = candidate_cost;
      closest_negative_mean = candidate_root;
    }
    if(NEWTON_STEPS <= ++step){
      return log((closest_positive_mean + closest_negative_mean)/2);
    }
    deriv = PoissonDeriv(candidate_root);
    possibly_outside = candidate_root - candidate_cost/deriv;
    if(possibly_outside < optimal_mean){
      if(closest_negative_cost==-INFINITY){
        throw 1;
      }
      return log((closest_positive_mean + closest_negative_mean)/2);
    }else{
      candidate_root = possibly_outside;
    }
  }while(NEWTON_EPSILON < ABS(candidate_cost));
  return log(candidate_root);
}

double PoissonLossPieceLog::get_smaller_root(double equals){
  double optimal_log_mean = argmin();
  double optimal_cost = getCost(optimal_log_mean);
  double left_cost = getCost(min_log_mean);
  if((equals < left_cost && left_cost < optimal_cost) ||
     (equals > left_cost && left_cost > optimal_cost)){
    return min_log_mean-1;
  }
  double candidate_root = optimal_log_mean - 1;
  double candidate_cost, possibly_outside, deriv;
  double closest_positive_cost = INFINITY, closest_positive_log_mean;
  double closest_negative_cost = -INFINITY, closest_negative_log_mean;
  if(optimal_cost < 0){
    closest_negative_cost = optimal_cost;
    closest_negative_log_mean = optimal_log_mean;
  }else{
    closest_positive_cost = optimal_cost;
    closest_positive_log_mean = optimal_log_mean;
  }
  int step=0;
  do{
    candidate_cost = getCost(candidate_root) - equals;
    if(0 < candidate_cost && candidate_cost < closest_positive_cost){
      closest_positive_cost = candidate_cost;
      closest_positive_log_mean = candidate_root;
    }
    if(closest_negative_cost < candidate_cost && candidate_cost < 0){
      closest_negative_cost = candidate_cost;
      closest_negative_log_mean = candidate_root;
    }
    if(NEWTON_STEPS <= ++step){
      return (closest_positive_log_mean + closest_negative_log_mean)/2;
    }
    deriv = getDeriv(candidate_root);
    possibly_outside = candidate_root - candidate_cost/deriv;
    if(possibly_outside < optimal_log_mean){
      candidate_root = possibly_outside;
    }else{
      return (closest_positive_log_mean + closest_negative_log_mean)/2;
    }
  }while(NEWTON_EPSILON < ABS(candidate_cost));
  return candidate_root;
}

void PoissonLossPieceLog::print(){
  Rprintf("%.10e %.10e %.10e %15f %15f %15f %d\n",
          Linear, Log, Constant,
          min_log_mean, max_log_mean, prev_log_mean, data_i);
}

// --- PiecewisePoissonLossLog methods (unchanged from reference) ---

void PiecewisePoissonLossLog::add(double Linear, double Log, double Constant){
  for(auto it=piece_list.begin(); it != piece_list.end(); it++){
    it->Linear += Linear;
    it->Log += Log;
    it->Constant += Constant;
  }
}

void PiecewisePoissonLossLog::multiply(double x){
  for(auto it=piece_list.begin(); it != piece_list.end(); it++){
    it->Linear *= x;
    it->Log *= x;
    it->Constant *= x;
  }
}

void PiecewisePoissonLossLog::set_prev_seg_end(int prev_seg_end){
  for(auto it=piece_list.begin(); it != piece_list.end(); it++){
    it->data_i = prev_seg_end;
  }
}

void PiecewisePoissonLossLog::findMean
(double log_mean, int *seg_end, double *prev_log_mean){
  for(auto it=piece_list.begin(); it != piece_list.end(); it++){
    if(it->min_log_mean <= log_mean && log_mean <= it->max_log_mean){
      *seg_end = it->data_i;
      *prev_log_mean = it->prev_log_mean;
      return;
    }
  }
  Rcpp::stop("findMean: log_mean=%e outside all piece intervals", log_mean);
}

double PiecewisePoissonLossLog::findCost(double mean){
  for(auto it=piece_list.begin(); it != piece_list.end(); it++){
    if(it->min_log_mean <= mean && mean <= it->max_log_mean){
      return it->getCost(mean);
    }
  }
  return INFINITY;
}

void PiecewisePoissonLossLog::Minimize
(double *best_cost, double *best_log_mean,
 int *data_i, double *prev_log_mean){
  double candidate_cost, candidate_log_mean;
  *best_cost = INFINITY;
  for(auto it=piece_list.begin(); it != piece_list.end(); it++){
    candidate_log_mean = it->argmin();
    if(candidate_log_mean < it->min_log_mean){
      candidate_log_mean = it->min_log_mean;
    }else if(it->max_log_mean < candidate_log_mean){
      candidate_log_mean = it->max_log_mean;
    }
    candidate_cost = it->getCost(candidate_log_mean);
    if(candidate_cost < *best_cost){
      *best_cost = candidate_cost;
      *best_log_mean = candidate_log_mean;
      *data_i = it->data_i;
      *prev_log_mean = it->prev_log_mean;
    }
  }
}

void PiecewisePoissonLossLog::print(){
  Rprintf("%10s %10s %15s %15s %15s %15s %s\n",
          "Linear", "Log", "Constant",
          "min_log_mean", "max_log_mean", "prev_log_mean", "data_i");
  for(auto it=piece_list.begin(); it != piece_list.end(); it++){
    it->print();
  }
}

// replaces set_to_min_less_of / set_to_min_more_of from the constrained solver.
// just call Minimize() and emit a flat line at the global min cost.
// (cf. funPieceListLog.cpp for the constrained versions)

void PiecewisePoissonLossLog::set_to_unconstrained_min_of
(PiecewisePoissonLossLog *input){
  double best_cost, best_log_mean, prev_log_mean_val;
  int prev_seg_end;
  input->Minimize(&best_cost, &best_log_mean, &prev_seg_end, &prev_log_mean_val);
  piece_list.clear();
  piece_list.emplace_back(
    0.0, 0.0, best_cost,
    input->piece_list.front().min_log_mean,
    input->piece_list.back().max_log_mean,
    prev_seg_end, best_log_mean);
}

// --- min-envelope helpers (from funPieceListLog.cpp, unchanged) ---

bool sameFuns
(PoissonLossPieceListLog::iterator it1,
 PoissonLossPieceListLog::iterator it2){
  return it1->Linear == it2->Linear &&
    it1->Log == it2->Log &&
    ABS(it1->Constant - it2->Constant) < NEWTON_EPSILON;
}

void PiecewisePoissonLossLog::push_piece
(PoissonLossPieceListLog::iterator it, double min_log_mean, double max_log_mean){
  if(max_log_mean <= min_log_mean) return;
  PoissonLossPieceListLog::iterator last=piece_list.end();
  --last;
  if(piece_list.size() &&
     sameFuns(last, it) &&
     it->prev_log_mean == last->prev_log_mean &&
     it->data_i == last->data_i){
    last->max_log_mean = max_log_mean;
  }else{
    piece_list.emplace_back(
      it->Linear, it->Log, it->Constant,
      min_log_mean, max_log_mean,
      it->data_i, it->prev_log_mean);
  }
}

void PiecewisePoissonLossLog::push_min_pieces
(PiecewisePoissonLossLog *fun1,
 PiecewisePoissonLossLog *fun2,
 PoissonLossPieceListLog::iterator it1,
 PoissonLossPieceListLog::iterator it2,
 int verbose){
  bool same_at_left;
  double last_min_log_mean;
  PoissonLossPieceListLog::iterator prev2 = it2;
  prev2--;
  PoissonLossPieceListLog::iterator prev1 = it1;
  prev1--;
  if(it1->min_log_mean < it2->min_log_mean){
    same_at_left = sameFuns(prev2, it1);
    last_min_log_mean = it2->min_log_mean;
  }else{
    last_min_log_mean = it1->min_log_mean;
    if(it2->min_log_mean < it1->min_log_mean){
      same_at_left = sameFuns(prev1, it2);
    }else{
      if(it1==fun1->piece_list.begin() &&
         it2==fun2->piece_list.begin()){
        same_at_left = false;
      }else{
        same_at_left = sameFuns(prev1, prev2);
      }
    }
  }
  PoissonLossPieceListLog::iterator next2 = it2;
  next2++;
  PoissonLossPieceListLog::iterator next1 = it1;
  next1++;
  bool same_at_right;
  double first_max_log_mean;
  if(it1->max_log_mean < it2->max_log_mean){
    same_at_right = sameFuns(next1, it2);
    first_max_log_mean = it1->max_log_mean;
  }else{
    first_max_log_mean = it2->max_log_mean;
    if(it2->max_log_mean < it1->max_log_mean){
      same_at_right = sameFuns(it1, next2);
    }else{
      if(next1==fun1->piece_list.end() &&
         next2==fun2->piece_list.end()){
        same_at_right = false;
      }else{
        same_at_right = sameFuns(next1, next2);
      }
    }
  }
  if(last_min_log_mean == first_max_log_mean) return;
  if(sameFuns(it1, it2)){
    push_piece(it1, last_min_log_mean, first_max_log_mean);
    return;
  }
  PoissonLossPieceLog diff_piece(
    it1->Linear - it2->Linear,
    it1->Log - it2->Log,
    it1->Constant - it2->Constant,
    last_min_log_mean, first_max_log_mean, -5, false);
  double mid_mean = (exp(first_max_log_mean) + exp(last_min_log_mean))/2;
  double cost_diff_mid = diff_piece.getCost(log(mid_mean));
  if(same_at_left && same_at_right){
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_log_mean, first_max_log_mean);
    }else{
      push_piece(it2, last_min_log_mean, first_max_log_mean);
    }
    return;
  }
  if(diff_piece.Log == 0){
    if(diff_piece.Linear == 0){
      if(diff_piece.Constant < 0){
        push_piece(it1, last_min_log_mean, first_max_log_mean);
      }else{
        push_piece(it2, last_min_log_mean, first_max_log_mean);
      }
      return;
    }
    if(diff_piece.Constant == 0){
      if(diff_piece.Linear < 0){
        push_piece(it1, last_min_log_mean, first_max_log_mean);
      }else{
        push_piece(it2, last_min_log_mean, first_max_log_mean);
      }
      return;
    }
    double log_mean_at_equal_cost = log(-diff_piece.Constant/diff_piece.Linear);
    if(last_min_log_mean < log_mean_at_equal_cost &&
       log_mean_at_equal_cost < first_max_log_mean){
      if(0 < diff_piece.Linear){
        push_piece(it1, last_min_log_mean, log_mean_at_equal_cost);
        push_piece(it2, log_mean_at_equal_cost, first_max_log_mean);
      }else{
        push_piece(it2, last_min_log_mean, log_mean_at_equal_cost);
        push_piece(it1, log_mean_at_equal_cost, first_max_log_mean);
      }
      return;
    }
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_log_mean, first_max_log_mean);
    }else{
      push_piece(it2, last_min_log_mean, first_max_log_mean);
    }
    return;
  }
  double cost_diff_left = diff_piece.getCost(last_min_log_mean);
  double cost_diff_right = diff_piece.getCost(first_max_log_mean);
  bool two_roots = diff_piece.has_two_roots(0.0);
  double smaller_log_mean, larger_log_mean;
  if(two_roots){
    smaller_log_mean = diff_piece.get_smaller_root(0.0);
    larger_log_mean = diff_piece.get_larger_root(0.0);
  }
  if(same_at_right){
    if(two_roots){
      double log_mean_at_crossing = smaller_log_mean;
      double log_mean_at_optimum = diff_piece.argmin();
      if(last_min_log_mean < log_mean_at_crossing &&
         log_mean_at_crossing < log_mean_at_optimum &&
         log_mean_at_optimum < first_max_log_mean){
        if(cost_diff_left < 0){
          push_piece(it1, last_min_log_mean, log_mean_at_crossing);
          push_piece(it2, log_mean_at_crossing, first_max_log_mean);
        }else{
          push_piece(it2, last_min_log_mean, log_mean_at_crossing);
          push_piece(it1, log_mean_at_crossing, first_max_log_mean);
        }
        return;
      }
    }
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_log_mean, first_max_log_mean);
    }else{
      push_piece(it2, last_min_log_mean, first_max_log_mean);
    }
    return;
  }
  if(same_at_left){
    if(two_roots){
      double log_mean_at_crossing = larger_log_mean;
      double log_mean_at_optimum = diff_piece.argmin();
      if(last_min_log_mean < log_mean_at_optimum &&
         log_mean_at_optimum < log_mean_at_crossing &&
         log_mean_at_crossing < first_max_log_mean){
        if(cost_diff_right < 0){
          push_piece(it2, last_min_log_mean, log_mean_at_crossing);
          push_piece(it1, log_mean_at_crossing, first_max_log_mean);
        }else{
          push_piece(it1, last_min_log_mean, log_mean_at_crossing);
          push_piece(it2, log_mean_at_crossing, first_max_log_mean);
        }
        return;
      }
    }
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_log_mean, first_max_log_mean);
    }else{
      push_piece(it2, last_min_log_mean, first_max_log_mean);
    }
    return;
  }
  double first_log_mean = INFINITY, second_log_mean = INFINITY;
  if(two_roots){
    bool larger_inside =
      last_min_log_mean < larger_log_mean && larger_log_mean < first_max_log_mean;
    bool smaller_inside =
      last_min_log_mean < smaller_log_mean &&
      0 < exp(smaller_log_mean) &&
      smaller_log_mean < first_max_log_mean;
    if(larger_inside){
      if(smaller_inside && smaller_log_mean < larger_log_mean){
        first_log_mean = smaller_log_mean;
        second_log_mean = larger_log_mean;
      }else{
        first_log_mean = larger_log_mean;
      }
    }else{
      if(smaller_inside){
        first_log_mean = smaller_log_mean;
      }
    }
  }
  if(second_log_mean != INFINITY){
    double before_mean = (exp(last_min_log_mean) + exp(first_log_mean))/2;
    double cost_diff_before = diff_piece.getCost(log(before_mean));
    if(cost_diff_before < 0){
      push_piece(it1, last_min_log_mean, first_log_mean);
      push_piece(it2, first_log_mean, second_log_mean);
      push_piece(it1, second_log_mean, first_max_log_mean);
    }else{
      push_piece(it2, last_min_log_mean, first_log_mean);
      push_piece(it1, first_log_mean, second_log_mean);
      push_piece(it2, second_log_mean, first_max_log_mean);
    }
  }else if(first_log_mean != INFINITY){
    double before_mean = (exp(last_min_log_mean) + exp(first_log_mean))/2;
    double cost_diff_before = diff_piece.getCost(log(before_mean));
    double after_mean = (first_max_log_mean + first_log_mean)/2;
    double cost_diff_after = diff_piece.getCost(after_mean);
    if(cost_diff_before < 0){
      if(cost_diff_after < 0){
        push_piece(it1, last_min_log_mean, first_max_log_mean);
      }else{
        push_piece(it1, last_min_log_mean, first_log_mean);
        push_piece(it2, first_log_mean, first_max_log_mean);
      }
    }else{
      if(cost_diff_after < 0){
        push_piece(it2, last_min_log_mean, first_log_mean);
        push_piece(it1, first_log_mean, first_max_log_mean);
      }else{
        push_piece(it2, last_min_log_mean, first_max_log_mean);
      }
    }
  }else{
    double cost_diff;
    if(first_max_log_mean == INFINITY){
      cost_diff = diff_piece.getCost(last_min_log_mean+1);
    }else{
      if(ABS(cost_diff_mid) < NEWTON_EPSILON){
        cost_diff = cost_diff_right;
      }else{
        cost_diff = cost_diff_mid;
      }
    }
    if(cost_diff < 0){
      push_piece(it1, last_min_log_mean, first_max_log_mean);
    }else{
      push_piece(it2, last_min_log_mean, first_max_log_mean);
    }
  }
}

void PiecewisePoissonLossLog::set_to_min_env_of
(PiecewisePoissonLossLog *fun1, PiecewisePoissonLossLog *fun2, int verbose){
  PoissonLossPieceListLog::iterator
    it1 = fun1->piece_list.begin(),
    it2 = fun2->piece_list.begin();
  piece_list.clear();
  while(it1 != fun1->piece_list.end() &&
        it2 != fun2->piece_list.end()){
    push_min_pieces(fun1, fun2, it1, it2, verbose);
    double last_max_log_mean = piece_list.back().max_log_mean;
    if(it1->max_log_mean == last_max_log_mean) it1++;
    if(it2->max_log_mean == last_max_log_mean) it2++;
  }
}

// --- 1-state DP loop (simplified from PeakSegFPOPLog.cpp) ---

int PoissonFPOPunconstrainedLog
(int *data_vec, double *weight_vec, int data_count,
 double penalty,
 double *cost_vec,
 int *end_vec,
 double *mean_vec,
 int *intervals_vec){
  double min_log_mean=INFINITY, max_log_mean=-INFINITY;
  for(int data_i=0; data_i<data_count; data_i++){
    if(data_vec[data_i] < 0) return ERROR_NEGATIVE_DATA;
    double log_data = log((double)data_vec[data_i]);
    if(log_data < min_log_mean) min_log_mean = log_data;
    if(max_log_mean < log_data) max_log_mean = log_data;
  }
  if(min_log_mean == max_log_mean){
    return ERROR_MIN_MAX_SAME;
  }

  // was N x 2 in PeakSegFPOPLog; just N x 1 here since no up/down states
  std::vector<PiecewisePoissonLossLog> cost_model_mat(data_count);
  PiecewisePoissonLossLog *cost, *cost_prev;
  PiecewisePoissonLossLog min_prev_cost;
  int verbose = 0;
  double cum_weight_i = 0.0, cum_weight_prev_i;

  for(int data_i=0; data_i<data_count; data_i++){
    cost = &cost_model_mat[data_i];
    cum_weight_i += weight_vec[data_i];
    if(data_i==0){
      // same init as PeakSegFPOPLog.cpp line 41
      cost->piece_list.emplace_back(
        1.0, -data_vec[0], 0.0,
        min_log_mean, max_log_mean, -1, INFINITY);
    }else{
      // flat line at global min -- replaces set_to_min_less_of
      min_prev_cost.set_to_unconstrained_min_of(cost_prev);
      min_prev_cost.set_prev_seg_end(data_i-1);
      min_prev_cost.add(0.0, 0.0, penalty/cum_weight_prev_i);
      // min-envelope: keep whichever is cheaper, switch or continue
      cost->set_to_min_env_of(&min_prev_cost, cost_prev, verbose);
      // weighted update (same formula as PeakSegFPOPLog lines 130-137)
      cost->multiply(cum_weight_prev_i);
      cost->add(
        weight_vec[data_i],
        -data_vec[data_i]*weight_vec[data_i],
        0.0);
      cost->multiply(1.0/cum_weight_i);
    }
    cum_weight_prev_i = cum_weight_i;
    cost_prev = cost;
  }

  // fill cost_vec and intervals; reuse last Minimize for backtracking
  double best_cost, best_log_mean, prev_log_mean;
  int prev_seg_end;
  for(int i=0; i<data_count; i++){
    cost = &cost_model_mat[i];
    intervals_vec[i] = cost->piece_list.size();
    cost->Minimize(cost_vec+i, &best_log_mean, &prev_seg_end, &prev_log_mean);
  }
  best_cost = cost_vec[data_count-1];

  // backtrack from last data point (cf. PeakSegFPOPLog lines 147-179)
  for(int i=0; i<data_count; i++){
    mean_vec[i] = INFINITY;
    end_vec[i] = -2;
  }
  mean_vec[0] = exp(best_log_mean);
  end_vec[0] = prev_seg_end;
  int out_i=1;
  while(0 <= prev_seg_end && out_i < data_count){
    cost = &cost_model_mat[prev_seg_end];
    if(prev_log_mean != INFINITY){
      best_log_mean = prev_log_mean;
    }
    cost->findMean(best_log_mean, &prev_seg_end, &prev_log_mean);
    mean_vec[out_i] = exp(best_log_mean);
    end_vec[out_i] = prev_seg_end;
    out_i++;
  }
  return 0;
}

// [[Rcpp::export]]
Rcpp::List PoissonFPOP(Rcpp::IntegerVector data_vec, double penalty){
  int n = data_vec.size();
  if(n < 2) Rcpp::stop("Need at least 2 data points");
  if(penalty < 0) Rcpp::stop("Penalty must be non-negative");

  std::vector<double> weight_vec(n, 1.0);
  std::vector<double> cost_vec(n);
  std::vector<int> end_vec(n);
  std::vector<double> mean_vec(n);
  std::vector<int> intervals_vec(n);

  int status = PoissonFPOPunconstrainedLog(
    data_vec.begin(), weight_vec.data(), n, penalty,
    cost_vec.data(), end_vec.data(), mean_vec.data(), intervals_vec.data());

  if(status == ERROR_NEGATIVE_DATA){
    Rcpp::stop("Negative values not allowed for Poisson loss");
  }
  if(status == ERROR_MIN_MAX_SAME){
    Rcpp::stop("All data values are identical -- need at least 2 distinct values");
  }

  // end_vec is filled in reverse order; count valid entries
  int n_segments = 0;
  for(int i=0; i<n; i++){
    if(end_vec[i] == -2) break;
    n_segments++;
  }

  // flip to chronological
  Rcpp::NumericVector means(n_segments);
  Rcpp::IntegerVector ends(n_segments);
  for(int i=0; i<n_segments; i++){
    means[i] = mean_vec[n_segments - 1 - i];
    ends[i] = end_vec[n_segments - 1 - i];
  }

  // 1-based breakpoints for R (ends[0..n_segments-2] + 1)
  Rcpp::IntegerVector breaks_vec;
  if(n_segments > 1){
    breaks_vec = Rcpp::IntegerVector(n_segments - 1);
    for(int i=0; i<n_segments-1; i++){
      breaks_vec[i] = ends[i] + 1;
    }
  }

  // undo the 1/cum_weight normalization
  double total_cost = cost_vec[n-1] * n;

  return Rcpp::List::create(
    Rcpp::Named("mean") = means,
    Rcpp::Named("end") = ends,
    Rcpp::Named("breaks") = breaks_vec,
    Rcpp::Named("n.segments") = n_segments,
    Rcpp::Named("cost") = total_cost,
    Rcpp::Named("cost.vec") = Rcpp::NumericVector(cost_vec.begin(), cost_vec.end()),
    Rcpp::Named("intervals") = Rcpp::IntegerVector(intervals_vec.begin(), intervals_vec.end())
  );
}
