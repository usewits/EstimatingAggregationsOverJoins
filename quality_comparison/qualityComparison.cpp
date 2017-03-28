#include <iostream>
#include <cmath>
#include <cstring>
#include <assert.h>
#include <algorithm>
#include <numeric>
#include <set>
#include <functional>
#include "sampleJoins.h"

#define MILLION 1000000

using namespace std;

//Generate weights (n elements), with a selected skew ratio and number of discrete values.
vector<double> get_distribution(int n, double skew, double ratio = 0.0, double n_discrete = 0.0) {
    if(ratio == 0.0) {
        //Either the ratio or n_discrete has to be defined 
        assert(n_discrete != 0);
        ratio = n_discrete;
    }
    vector<double> w(n);
    for(int i=0; i<n; i++)
        w[i] = pow(mtwist_drand(mt), skew); //the weights are in [0,1[ with a (polynomial) skew
    double max_w = *max_element(w.begin(), w.end());
    double sum_w = 0;
    for(int i=0; i<n; i++) {
        w[i] = (w[i]/max_w)*(ratio-1.0)+1.0; //the weights are in [1, ratio[
        sum_w += w[i];
    }
    //normalise the weights
    for(int i=0; i<n; i++)
        w[i] /= sum_w;

    if(n_discrete > 0) {
        max_w = *max_element(w.begin(), w.end());
        for(int i=0; i<n; i++) {
            w[i] *= n_discrete/max_w;
            w[i] = round(w[i]);
        }
    }
    return w;    
}

void getValues(double& a, double& b, double& c, const tdd& t) {
    a = get<0>(t);
    b = get<1>(t);
    c = get<2>(t);
}


//O(1) time
double HWS_heuristic_simple(const vector<double>& w, double sigma, double k_factor, int m) {
    return (double)m*(double)m;
}

double HWS_heuristic_complete(const vector<double>& w, double w_ratio, double sigma, double k_factor, int m) {
    double w_max = (*max_element(w.begin(), w.end()));//O(|w|) time
    double w_min = (*min_element(w.begin(), w.end()));
    double sigma_factor = 1.0/log(1.0/sigma);
    return k_factor*sigma_factor* (double)m*(double)m * w_max/w_min;
}




//Generic function to estimate aggregates over joins
//It can be used to obtain SSJ, HSSJ, WS-Join, HWS-Join or US-Join estimates (both filtered and unfiltered)
//Runs in O(k+n2) time (not as fast as possible, in favour of shorter code)
//recompute_normalisation should be set to true whenever any of the following variables have changed: h1, h2, R1, R2, R1_filter, R2_filter
//    This adds O(n1) to the runtime
//recompute_cdf causes memoisation of the cdf on R1. This cdf is invalidated if the normalisation is recomputed.
//    This adds O(n1) to the runtime.
//    If no memoized cdf is available, it will be computed by the range_sampler if necessary instead, taking between O(1) and O(n1) time
    
//Uses O(n1) = ~3*n1*(2*64) bits of memory
//Output probability is h1*h2
double generic_sample_join(function<double(double,double)> h1, function<double(double)> h2, int m,
                                const vector<pdd>& R1, const vector<pdd>& R2,
                                function<vector<int>(int, const vector<double>&, vector<double>*)> range_sampler,
                                            //sample(sample_size, weights)
                                function<double(double, double, double)> aggregation_f,
                                function<bool(double, double)> R1_filter,//Ri_filter are predicates; true => selected
                                function<bool(double, double)> R2_filter,
                                bool filtered_estimator, double filter_selectivity,
                                bool recompute_normalisation, bool recompute_cdf) {

    //Compute (filtered) stratum weights and cdfs (O(n2) time, O(n2) memory)
    Tstrat R2_stratified = stratify(R2);//O(n2) memory
    map<double, double> R2_stratum_weights;
    map<double, double> R2_filtered_stratum_weights;
    map<double, vector<double> > R2_stratum_cdfs;//O(n2) memory
    for(auto stratum : R2_stratified) {//O(n2) time
        double key = stratum.first;
        double norm = 0.0;
        double filtered_norm = 0.0;
        vector<double> stratum_weights(stratum.second.size());
        for(int i=0; i<stratum.second.size(); i++) {
            pdd t2 = stratum.second[i];
            stratum_weights[i] = h2(t2.second);
            norm += stratum_weights[i];
            if(R2_filter(t2.first, t2.second))
                filtered_norm += stratum_weights[i]; 
        }
        R2_stratum_weights[key] = norm;
        R2_filtered_stratum_weights[key] = filtered_norm;
        R2_stratum_cdfs[key] = get_cdf(stratum_weights);
    }

    
    static double normalisation = 0.0;          //Total weight of all elements in J
    static double filtered_normalisation = 0.0; //Total weight of selection sigma(J)
    static vector<double> R1_sample_weights(R1.size());                //Sampling weights in R1 (n1 memory)
    static vector<double> R1_filtered_sample_weights(R1.size(), 0.0);  //Filtered sampling weights (n1 memory)
    static vector<double> *R1_sample_weights_cdf = NULL; //At first, no cdf is available

    if(recompute_normalisation || recompute_cdf) {
        if(R1_sample_weights_cdf != NULL) {//deallocate cdf if necessary
            vector<double> *tmp = R1_sample_weights_cdf;
            R1_sample_weights_cdf = NULL;
            delete tmp;
        }
    }
    if(recompute_normalisation) {
        //Compute normalisation factors (O(n1) time, 2*n1 memory)
        //These depend on: h1, h2, R1_filter, R2_filter, R1, R2 (and none of the other arguments)
        normalisation = 0.0;
        filtered_normalisation = 0.0;
        R1_sample_weights = vector<double>(R1.size());
        R1_filtered_sample_weights = vector<double>(R1.size(), 0.0);
        
        for(int i=0; i<R1.size(); i++) {//O(n1) time
            pdd t1 = R1[i];
            R1_sample_weights[i] = h1(t1.first, t1.second) * R2_stratum_weights[t1.first];
            normalisation += R1_sample_weights[i];
            if(R1_filter(t1.first, t1.second)) {
                R1_filtered_sample_weights[i] = h1(t1.first, t1.second) * R2_filtered_stratum_weights[t1.first];
                filtered_normalisation += R1_filtered_sample_weights[i];
            }
        }
        R1_sample_weights_cdf = NULL; //invalidate cdf (it depends on the normalisation)
    }

    if(recompute_cdf) {
        R1_sample_weights_cdf = new vector<double>(get_cdf(R1_sample_weights));
    }
    
    //Construct sample (O(k+m'[+n1]) time, O(k) memory)
    double over_sampling_factor = 1.2;
    int over_sampling_constant = 100;
    int S_size = round(over_sampling_constant+ceil(over_sampling_factor*m/filter_selectivity));
    vector<int> S_indices = range_sampler(S_size, R1_sample_weights, R1_sample_weights_cdf);
                            //full HWS heuristic: O(n1) time, O(k) memory
                            //simple HWS heuristic: O(k) time and memory
                            //Reason: min and max of R1_sample_weights are not memoized
    vector<pdd> S(S_size);
    vector<double> S_weights(S_size);
    for(int i=0; i<S_size; i++) {//O(m'=m/selectivity)=O(S_size) time
        S[i] = R1[S_indices[i]];
        S_weights[i] = R1_sample_weights[S_indices[i]];
    }
    vector<tdd> sample = minijoin(S, R2_stratified);//O(m') time and memory
    
    int filtered_sample_size = 0;

    for(auto t : sample) {//O(m') time
        double tA, tB, tC;
        getValues(tA, tB, tC, t);
        if(R1_filter(tA, tB) && R2_filter(tA, tC)) {
            filtered_sample_size++;
        }
    }

    //Reduce sample size until the filtered_sample_size equals m (O(m') time)
    assert(filtered_sample_size >= m);//If this is not the case, S_size is too small
    while(filtered_sample_size > m) {
        double tA, tB, tC;
        getValues(tA, tB, tC, sample.back());
        if(R1_filter(tA, tB) && R2_filter(tA, tC)) {
            filtered_sample_size--;
        }
        sample.pop_back();
    }
    
    //Compute estimate (O(m') time)
    double estimate = 0.0;
    for(auto t : sample) {//O(m') time
        double tA, tB, tC;
        getValues(tA, tB, tC, t);
        double w_t = h1(tA, tB)*h2(tC);//non normalised weights
        if(R1_filter(tA, tB) && R2_filter(tA, tC)) {
            estimate += aggregation_f(tA, tB, tC)/w_t;
        }
    }
    
    if(filtered_estimator) {//Correct for filter using filter-specific normalisation
        estimate *= filtered_normalisation/(double) filtered_sample_size;
    } else {//Use default normalisation (if filter is used, convergence is not guaranteed)
        estimate *= normalisation/(double) sample.size();
    }
    return estimate;
}


//total memory requirement: ~ 11*n1*64 bits
int main() {
    //initialize rng
    mt = mtwist_new();
    //mtwist_seed(mt, 832982837UL);
    mtwist_seed(mt, time(NULL));

    int m = 100;
    double k_factor = 1.0;
    double sigma = 0.99;
    double sigma_factor = 1.0/log(1/sigma);

    // Generate R1
    int n1          = 200*MILLION;
    double skew1    = 1.0;
    double ratio1   = 20.0;
    int n_discrete1 = 10.0;
    vector<pdd> R1;//~n1*(2*64) bits of memory
    {
        vector<double> R1A = get_distribution(n1,skew1,ratio1,n_discrete1);
        vector<double> R1B = get_distribution(n1,1.0,n1);
        R1 = zipvec(R1A, R1B);
    }
    auto stratR1 = stratify(R1);//~n1*(2*64) bits of memory

    // Generate R2
    int n2          = 2000;
    double skew2    = 1.0;
    double ratio2   = 50.0;
    int n_discrete2 = 10.0;
    vector<pdd> R2;
    {
        vector<double> R2A = get_distribution(n2,skew2,ratio2,n_discrete2);
        vector<double> R2C = get_distribution(n2,1.0,n2);
        R2 = zipvec(R2A, R2C);
    }
    Tstrat stratR2 = stratify(R2);
   
    // Define lambda functions
    auto aggregate_f = [] (double A, double B, double C) -> double {return C;};

    auto h1_unif =     []   (double A, double B) -> double {return 1.0;};
    auto h1_US = [&stratR2] (double A, double B) -> double {return 1.0/(double)(stratR2[A].size());};
    auto h1_weighted = []   (double A, double B) -> double {return 1.0;};

    auto h2_unif = [] (double C) -> double {return 1.0;};
    auto h2_weighted = [] (double C) -> double {return C;};

    auto HWS_heuristic = *HWS_heuristic_simple;

    auto exact_sampler = [] (int m, const vector<double>& w, vector<double>* c_w) -> vector<int> {
                                bool recompute_c_w = (c_w == NULL);
                                if(recompute_c_w)
                                    c_w = new vector<double>(get_cdf(w));//O(|w|) time
                                vector<int> result = weighted_sample_indices(w.size(), *c_w, m);
                                if(recompute_c_w)
                                    delete c_w;
                                return result;
                            };
    auto heuristic_sampler = [&sigma,&k_factor,&HWS_heuristic] (int m, const vector<double>& w, vector<double>* c_w) -> vector<int> {

                                int k = round(HWS_heuristic(w, sigma, k_factor, m));//O(1) or O(|w|) time
    
                                vector<int> U = sample_indices(w.size(), k);//O(k) time
                                vector<double> U_w(k);
                                for(int i=0; i<k; i++) {//O(k) time
                                    U_w[i] = w[U[i]];
                                }
                                return weighted_sample(U, get_cdf(U_w), m);//O(k) time
                            };

    //Selection filters: tuples that produce a true are selected
    auto no_filter = [] (double X, double Y) -> bool {return true;};
    auto range_filter = [] (double X, double Y) -> bool {return Y > 1.5e-3;};
    auto rand_filter = [] (double X, double Y) -> bool {int* Z; Z=(int*)&Y; return (*Z)%2;};

    auto R1_filter = no_filter;
    auto R2_filter = rand_filter;
 
    
    //Compute the sampling weights required for SSJ
    vector<double> SSJ_prob(n1);//~n1*64 bits of memory
    for(int i=0; i<n1; i++) {
        double key = R1[i].first;
        SSJ_prob[i] = stratR2[key].size();
           //note stratR2[key].size() = m_2(t_1.A)
    }

    //A list of generic-sample-join parameters and their names
    set<int> sampling_methods_used = {0,1,2,3,4};
    string                          sample_types[] = {"SSJ     ","HSSJ    ","WS-Join ","HWS-Join",  "US-Join "};
    function<double(double,double)> h1_functions[] = { h1_unif,   h1_unif,  h1_weighted,h1_weighted, h1_US};
    function<double(double)>        h2_functions[] = { h2_unif,   h2_unif,  h2_weighted,h2_weighted, h2_unif};
    bool                            is_heuristic[] = {   false,      true,        false,       true,   false};
    function<vector<int>(int, const vector<double>&, vector<double>*)> samplers[] = 
                        { exact_sampler, heuristic_sampler, exact_sampler, heuristic_sampler, exact_sampler};
   
    //A list of filter methods and their names
    set<int> filter_methods_used = {0,1,2};
    string filter_types[] = {"full      ","filtered  ","fltr.naive"};
    function<bool(double,double)> R1_filters[] = {no_filter, R1_filter, R1_filter};
    function<bool(double,double)> R2_filters[] = {no_filter, R2_filter, R2_filter};
    bool filtered_estimations[] = {false, true, false};

    //Compute and print the true aggregate values for each filter mode (actually the same for filtered and fltr.naive)

    vector<double> true_aggregates(3, 0.0);
    vector<long long> filtered_join_size(3, 0);
    vector<double> selectivities(3, 0.0);

    bool aggregate_f_independent_of_B = true;//True aggregate can be computed faster if simple.

    long long full_join_size = 0;
    if(aggregate_f_independent_of_B) {//O(n1+n2) time
        map<double, double> R2_exact_aggregates[3];
        map<double, int> R2_exact_sizes[3];
        for(auto strat2 : stratR2) {
            double tB = -9999;
            for(auto t2 : strat2.second) {
                double tA = t2.first;
                double tC = t2.second;
                for(int i_f : filter_methods_used) {
                    if(R1_filters[i_f](tA, tB) && R2_filters[i_f](tA, tC)) {
                        R2_exact_aggregates[i_f][tA] += aggregate_f(tA, tB, tC);
                        R2_exact_sizes[i_f][tA]++;
                    }
                }
            }
        }
        
        for(auto strat1 : stratR1) {
            double a = strat1.first;
            auto strat2it = stratR2.find(a);
            if(strat2it == stratR2.end())
                continue; //key does not join
            for(auto t1 : strat1.second) {
                double tA = t1.first;
                double tB = t1.second;
                full_join_size+= R2_exact_sizes[0][tA];

                for(int i_f : filter_methods_used) {
                    if(R1_filters[i_f](tA, tB)) {
                        true_aggregates[i_f] += R2_exact_aggregates[i_f][tA];
                        filtered_join_size[i_f] += R2_exact_sizes[i_f][tA];
                    }
                }
            }
        }
    } else {//O(|J|) time
        for(auto strat1 : stratR1) {
            double a = strat1.first;
            auto strat2it = stratR2.find(a);
            if(strat2it == stratR2.end())
                continue; //key does not join
            for(auto t1 : strat1.second)
            for(auto t2 : strat2it->second) {
                tdd j = make_tuple(t1.first, t1.second, t2.second);//Here j : J where J the full join
                                                                   //J = join(stratR1, stratR2);
                full_join_size++;

                double tA, tB, tC;
                getValues(tA, tB, tC, j);

                for(int i_f : filter_methods_used) {
                    if(R1_filters[i_f](tA, tB) && R2_filters[i_f](tA, tC)) {
                        true_aggregates[i_f] += aggregate_f(tA, tB, tC);
                        filtered_join_size[i_f] ++;
                    }
                }
            }
        }
    }
    cout << "Join size: " << full_join_size << endl;

    for(int i_f : filter_methods_used) {
        selectivities[i_f] = filtered_join_size[i_f]/(double) full_join_size;
        cout << "Exact aggregation (" << filter_types[i_f] << ") :" << true_aggregates[i_f] << " (selectivity " << selectivities[i_f]*100 << "% -> sample size ~ "<< round(m/selectivities[i_f])<< ")" << endl;
    }
    
     
    //THE EXPERIMENTS

    int nruns = 1000;
    while(true) {
        cout << "Sample size? m = " << flush;
        cin >> m;
        cout << "Runs? nruns = " << flush;
        cin >> nruns;
        //Couls also set HWS mode..
        
        map< pair<int, int>, vector<double> > relative_errors;
        
        //Initialize the relative_errors object
        for(int i_s : sampling_methods_used)
        for(int i_f : filter_methods_used) {
            relative_errors[make_pair(i_s, i_f)] = vector<double>(nruns, 0);
        }
 
        cout << "Running " << nruns << "*" << relative_errors.size() << " experiments..." << endl;
 
        //Remove HWS-based methods if HWS causes oversampling 
        double k_dbl = HWS_heuristic(SSJ_prob, sigma, k_factor, m);
        cout << "k = " << k_dbl << " (should be smaller than " << R1.size() << " for AWS)"<< endl;
        if(k_dbl > R1.size()) {
            cout << "WARNING: Skipping Heuristic methods!" << endl;
            sampling_methods_used.erase(1);
            sampling_methods_used.erase(3);
        } else {
            sampling_methods_used.insert(1);
            sampling_methods_used.insert(3);
        }

        
        char confirm = 'n';
        cout << "Are you sure? (y/n) " << flush;
        cin >> confirm;
        if(confirm != 'y')
            continue;
            
        //Run experiments for each setting nruns times
        
        for(int i_f : filter_methods_used)
        for(int i_s : sampling_methods_used) {
            int progress_width = 50;
            for(int run_i=0; run_i<nruns; run_i++) {
                if(floor(progress_width*(run_i+1)/(double)nruns) > floor(progress_width*(run_i)/(double)nruns)) {
                    int n_bars = round(progress_width*run_i/(double)nruns);//out of 100
                    cout << " [";
                    for(int progress = 0; progress < progress_width; progress++) {
                        if(progress < n_bars)
                            cout << "#";
                        else
                            cout << " ";
                    }
                    if(progress_width == n_bars)
                        cout << "] DONE! " << endl;
                    else
                        cout << "] " << round(100*run_i/(double)nruns) << "%\r" << flush;
                }
                bool recompute_normalisation = (run_i == 0);
                bool recompute_cdf = recompute_normalisation && !is_heuristic[i_s];
                    //Make generic_sample_join memoise normalisation only if it is not a heuristic sample join
                    //since heuristic sample joins do not require the full cdf
                double estimate = generic_sample_join(h1_functions[i_s], h2_functions[i_s], 
                                                      m, R1, R2, samplers[i_s], aggregate_f, 
                                                      R1_filters[i_f], R2_filters[i_f], filtered_estimations[i_f],
                                                      selectivities[i_f], recompute_normalisation, recompute_cdf);
                relative_errors[make_pair(i_s, i_f)][run_i] = abs(true_aggregates[i_f]-estimate)/true_aggregates[i_f];
            }

            //Print the results (CI intervals)
            cout << sample_types[i_s] << "(" << filter_types[i_f] << "):" << endl;
            //show_sigma_levels(relative_errors[make_pair(i_s, i_f)]);
            show_sd(relative_errors[make_pair(i_s, i_f)]);

        }
    }

    return 0;
}
