#include <Rcpp.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector DAMOCLES_sim_cpp(Rcpp::DataFrame pa, double gamma, double mu) {
  
  double tstep = 0, wt, event;
  IntegerVector extant = pa["extant"];
  IntegerVector state = pa["state"];
  std::vector<double> start = pa["age.start"];
  std::vector<double> end = pa["age.end"];
  std::vector<int> parent = pa["p"];
  std::vector<int> daugther = pa["d"];
  
  double tlastevent = *std::max_element(end.begin(), end.end());
  
  int nExt = sum(extant);
  int nP = sum(state);
  int nA = nExt - nP;
  std::vector<int> focP, focA;
  
  double tnextevent = tlastevent;
  for(auto j = 0; j < extant.size(); j++){
    if(extant[j] == 1){
      if(state[j] == 1){
        focP.push_back(j);
      }else{
        focA.push_back(j);
      }
      if(end[j] < tnextevent){
        tnextevent = end[j];
      }
    }
  }
  // rate totals across extant species
  double gamma_tot = gamma * (nExt - nP);
  double mu_tot = mu * nP;
  // maximum rates based on the number of extant species
  double mu_max = mu * nExt;
  double gamma_max = gamma * nExt;
  double tot_max = gamma_max + mu_max;
  // relative rates based on the ratio between the actual rates and the maximum rates
  double mu_rel = mu_tot / tot_max;
  double gamma_rel = gamma_tot / tot_max;
  
  // create a vector with the probabilities of each event, equal to their total rate
  NumericVector prob = NumericVector::create(1 - (gamma_rel + mu_rel), gamma_rel, mu_rel);
  
  if(tot_max <= 0){
    // if the total rate is <= 0 no event will be possible so kill the simulation
    Rcout << "tot rate error" << std::endl;
    return -1;
  }
  
  while(tstep < tlastevent){
    wt = Rcpp::rexp(1, tot_max)[0];
    tstep += wt;
    
    if(tstep < tnextevent){
      // sample an event
      event = Rcpp::sample(3, 1, false, prob)[0] - 1;
      
      // colonisation
      if(event == 1){
        int foc = Rcpp::sample(nA, 1, false)[0] - 1;
        focP.push_back(focA[foc]);
        std::vector<int>::iterator nth = focA.begin() + foc;
        focA.erase(nth);
        nP++;
        nA--;
        mu_tot = mu * nP;
        gamma_tot = gamma * nA;
        mu_rel = mu_tot / tot_max;
        gamma_rel = gamma_tot / tot_max;
        prob[0] = 1 - (gamma_rel + mu_rel);
        prob[1] = gamma_rel;
        prob[2] = mu_rel;
      }
      // local extinction
      if(event == 2){
        int foc = Rcpp::sample(nP, 1, false)[0] - 1;
        focA.push_back(focP[foc]);
        std::vector<int>::iterator nth = focP.begin() + foc;
        focP.erase(nth);
        nP--;
        nA++;
        mu_tot = mu * nP;
        gamma_tot = gamma * nA;
        mu_rel = mu_tot / tot_max;
        gamma_rel = gamma_tot / tot_max;
        prob[0] = 1 - (gamma_rel + mu_rel);
        prob[1] = gamma_rel;
        prob[2] = mu_rel;
      }
      // if it's time for the next event, update the states, presences/absences, and focal species
    }else if (tnextevent != tlastevent){
      // find which species are speciating (might me more than one so a while loop is needed)
      std::vector<double>::iterator it = end.begin();
      while ((it = std::find_if(it, end.end(), [&tnextevent](double x){return x == tnextevent; })) != end.end()){
        int foc = std::distance(end.begin(), it);
        // set the species to extinct
        extant[foc] = 0;
        std::vector<int>::iterator iter1 = std::find(focA.begin(), focA.end(), foc);
        std::vector<int>::iterator iter2 = std::find(focP.begin(), focP.end(), foc);
        // if absent, change state to absent as well
        if (iter1 != focA.end()){ // == if iter1 != focA.end() this means that it is absent
          state[foc] = 0;
          focA.erase(iter1);
        }else{
          state[foc] = 1;
          focP.erase(iter2);
        }
        // find if the species has any daughters
        if(std::find(parent.begin(), parent.end(), daugther[foc]) != parent.end()){
          std::vector<int> B;
          int da = daugther[foc]; // name of the ancestral species
          // we're loopin through the parents of the daughters to see if any has the ancestral species as its parent
          std::vector<int>::iterator pa = parent.begin();
          while ((pa = std::find_if(pa, parent.end(), [&da](int x){return x == da; })) != parent.end()){
            int dau = std::distance(parent.begin(), pa);
            extant[dau] = 1;
            B.push_back(dau);
            pa++;
          }
          // set the states for the daugther lineages, which depend on the state of the parent (either 1 or 0 if in the comm or both 0 if p is absent)
          if(state[foc] == 1){
            IntegerVector A = {0,1};
            IntegerVector sts = Rcpp::sample(A, 2, false);
            
            int z = 0;
            for(std::vector<int>::const_iterator k = B.begin(); k != B.end(); ++k) {
              if(sts[z] == 1){
                focP.push_back(*k);
              }else{
                focA.push_back(*k);
              }
              z++;
            }
            // if the parent is absent, both daughters will be as well
          }else{
            for(std::vector<int>::const_iterator k = B.begin(); k != B.end(); ++k) {
              focA.push_back(*k);
            }
          }
        }
        it++;
      }
      // update the counters
      nExt = sum(extant);
      nP = focP.size();
      nA = focA.size();
      mu_tot = mu * nP;
      gamma_tot = gamma * nA;
      
      mu_max = mu * nExt;
      gamma_max = gamma * nExt;
      tot_max = gamma_max + mu_max;
      
      mu_rel = mu_tot / tot_max;
      gamma_rel = gamma_tot / tot_max;
      prob = {1 - (gamma_rel + mu_rel), gamma_rel, mu_rel};
      // change the tstep to the time of speciation
      tstep = tnextevent;
      // determine when the next event is going to happen
      tnextevent = tlastevent;
      for(int j = 0; j < extant.size(); j++){
        if(extant[j] == 1){
          if(end[j] < tnextevent){
            tnextevent = end[j];
          }
        }
      }
    }
  }
  for(std::vector<int>::const_iterator l = focP.begin(); l != focP.end(); ++l) {
    state[*l] = 1;
  }
  return state;
}

