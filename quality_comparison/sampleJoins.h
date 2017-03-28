#include <cmath>
#include <map>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <stdlib.h>
#include <tuple>

#include "mtwist.h"

#define pdd pair<double, double>
#define tdd tuple<double, double, double>
#define Tstrat map<double, vector<pdd> >
using namespace std;

//The only global variable
mtwist* mt;

//input:  a weight vector (need not be normalised)
//output: corresponding CDF (from w[0]/W to 1)
vector<double> get_cdf(const vector<double>& w) {
    vector<double> result(w.size());
    result[0]=w[0];
    for(int i=1; i<w.size(); i++) //set result to cumsum of w
        result[i] = result[i-1]+w[i];
    for(int i=0; i<w.size(); i++) //normalize result (last element should be 1)
        result[i]/=result[w.size()-1];
    return result;
}

//input:  two columns of data
//output: data stratified by first column
Tstrat stratify(const vector<pdd>& R) {
    Tstrat result;
    for(vector<pdd>::const_iterator Rit = R.begin(); Rit!=R.end(); ++Rit) {
        if(result.find(Rit->first) == result.end()) {
            result.insert(make_pair(Rit->first,vector<pdd>()));//for first element in strata, insert empty vector
        }
        result[Rit->first].push_back(*Rit);
    }
    return result;
}

//input:  two columns of data
//output: map from first column to size of strata
map<double, int> stratify_counts(const vector<pdd>& R) {
    Tstrat Rstrat = stratify(R);
    
    map<double, int> result;
    for(Tstrat::const_iterator Rit = Rstrat.begin(); Rit!=Rstrat.end(); ++Rit) {
        result[Rit->first]=(Rit->second).size();
    }
    return result;
}

//Multiply the counts of strata (can be used to compute full join size)
map<double, int> strat_count_multiply(const map<double, int>& a, const map<double, int>& b) {
    map<double, int> result;
    for(map<double, int>::const_iterator ait = a.begin(); ait != a.end(); ++ait) {
        map<double, int>::const_iterator bit = b.find(ait->first);
        if(bit != b.end()) {
            result.insert(make_pair(ait->first, ait->second*bit->second));
        }
    }
    return result;
}

//Obtain sample with replacement of size k
template <typename T> vector<T> sample(const vector<T>& R, int k) {
    vector<T> result(k);
    int n = R.size();
    for(int i=0; i<k; i++) {
        result[i]=R[mtwist_uniform_int(mt,0,n-1)];
    }
    return result;
}

//Obtain sample with replacement of size k from a list of indices {0,1,...,n-2,n-1}
vector<int> sample_indices(int n, int k) {
    vector<int> result(k);
    for(int i=0; i<k; i++) {
        result[i]=mtwist_uniform_int(mt,0,n-1);
    }
    return result;
}

//input:   R   - population to sample from
//         c_p - normalized cumilative probability distribution 
//         m   - sample size
//output:  weighted sample with replacement of size m
template <typename T> vector<T> weighted_sample(const vector<T>& R, const vector<double>& c_p, int m) {
    vector<T> result(m);
    int n = R.size();
    for(int i=0; i<m; i++) {
        double random_variate = mtwist_drand(mt);
        vector<double>::const_iterator c_p_it = upper_bound(c_p.begin(), c_p.end(), random_variate);
        int index = c_p_it-c_p.begin();
        result[i]=R[index];
    }
    return result;
}

//same as weighted_sample, but R is replaced by {0, 1, ..., n-2, n-1}
vector<int> weighted_sample_indices(int n, const vector<double>& c_p, int m) {
    vector<int> result(m);
    for(int i=0; i<m; i++) {
        double random_variate = mtwist_drand(mt);
        vector<double>::const_iterator c_p_it = upper_bound(c_p.begin(), c_p.end(), random_variate);
        int index = c_p_it-c_p.begin();
        result[i]=index;
    }
    return result;
}

//same as weighted_sample, but heuristic.
//w_ratio = w_max / w_min
//delta is the desired maximum normalised absolute difference between the projected total weight from U and the total weight in R
template <typename T> vector<T> approximate_weighted_sample(const vector<T>& R, const vector<double>& c_p,
                                                            int m, double w_ratio, double delta) {
    double min_k = ceil(w_ratio*(double)m*(double)m);

    if(min_k > R.size()) {//NOTE: in practice we may want to switch from AWS to WS at much lower sampling fraction, depending on runtime cutoff.
                          //This switch is just here to avoid run-away runtimes, and will not influence determination of crossover since it is much beyond the crossover point
        cout << ("WARNING: oversampling in AWS! Early switch to regular sampling.\n");
        cout << ("AWS effective sampling fraction 1.0\n");
        return(weighted_sample(R, c_p, m));
    }

    double memory_factor = 2.0;
    double U_weight = 0;

    vector<double> w_U;
    vector<int> U_indices;
    w_U.reserve(ceil(memory_factor*min_k));
    U_indices.reserve(ceil(memory_factor*min_k));

    int n=R.size();
    int k;
    for(k=0; k<min_k || abs(U_weight*n/(double)k-1.0) > delta; k++) {
                        //Note; U_weight*n/(double)k is the projected total weight from U
        //int j = *sample_indices(n,1).begin();
        int j=mtwist_uniform_int(mt,0,n-1);
        U_indices.push_back(j);

        if(j == 0)
            w_U.push_back(c_p[j]);
        else
            w_U.push_back(c_p[j] - c_p[j-1]);

        U_weight += w_U[k];
        if(k > R.size()) {
            cout << ("WARNING: oversampling in AWS! Late switch to regular sampling.\n");
            cout << ("AWS effective sampling fraction 1.0\n");
            return(weighted_sample(R, c_p, m));
        }
    }
    if(k > min_k) {//Only print effective sampling fraction if it is not min_k/n1
        cout << "AWS effective sampling fraction " << k/(double)R.size() << " (default sampling fraction is " << min_k/(double)R.size() << ")" << endl;
    }

    vector<double> c_p_U = get_cdf(w_U);
    
    vector<int> S_indices = weighted_sample(U_indices, c_p_U, m);
    vector<T> S(m);
    for(int i=0; i<m; i++) {
        S[i] = R[S_indices[i]];
    }
    return S;
}


vector<pdd> zipvec(const vector<double>& l, const vector<double>& r) {
    assert(l.size() == r.size());
    int n=l.size();
    vector<pdd> result(n);
    for(int i=0; i<n; i++) {
        result[i] = make_pair(l[i],r[i]);
    }
    return result;
}

//Join two 2-column relations on the first attribute to form a 3-column relation
vector<tdd> join(const Tstrat& R1, const Tstrat& R2) {
    vector<tdd> result;
    for(auto strat1 : R1) {
        double a = strat1.first;
        auto strat2it = R2.find(a);
        if(strat2it == R2.end())
            continue; //key does not join
        for(auto t1 : strat1.second)
        for(auto t2 : strat2it->second) {
            result.push_back(make_tuple(t1.first, t1.second, t2.second));
        }
    }
    return result;
}

vector<pdd> printvec(const vector<pdd>& V) {
    for(auto v : V) {
        cout << v.first << "\t" << v.second << endl;
    }
    cout << endl;
}


vector<tdd> minijoin(const vector<pdd>& S, const Tstrat& R2) {
    vector<tdd> result;
    for(auto t1 : S) {
        auto strat2it = R2.find(t1.first);
        if(strat2it == R2.end())
            continue; //key does not join
        auto t2 = sample(strat2it->second,1)[0];
        result.push_back(make_tuple(t1.first, t1.second, t2.second));
    }
    return result;
}


void show_sigma_levels(vector<double>& relative_errors) {
    sort(relative_errors.begin(), relative_errors.end());
    double sigmas[] = {0.9, 0.95, 0.99};
    for(double sigma : sigmas) {
        double epsilon = relative_errors[(int)round(sigma*relative_errors.size())];
        cout << "\tapproximation (sigma = "<<sigma<<", epsilon = "<<epsilon*100<<"%)"<<endl;
    }
}

void show_sd(vector<double>& relative_errors) {
    double mean = accumulate(relative_errors.begin(),relative_errors.end(),0.0)/relative_errors.size();
    double variance = 0;
    for(auto err : relative_errors)
        variance += (mean-err)*(mean-err);
    double sd = sqrt(variance);
    cout << "\tapproximation (mean = "<<mean<<", sd = "<<sd<<")"<<endl;
}
