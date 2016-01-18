#ifndef __UTIL_HPP__
#define __UTIL_HPP__

#include <vector>
#include <stack>
#include <sstream> 

namespace util {
  using namespace std;

  // Note median actually modifies the order of the elements in the
  // given vector!! This happens because nth_element's O(N) selection
  // algorithm will swap elements around.
  template< class T > 
  struct _median {
    typedef typename vector<T>::iterator it_t;
    T operator () (it_t begin, it_t end) {
      it_t it = begin;
      size_t N = end - begin;
      if(N == 1) return *begin;
      if(N % 2 == 0) { // even number of elements
	it_t mid = begin + ( end - begin ) / 2 - 1;
	nth_element(begin, mid, end);
	it_t mid1 = begin + ( end - begin ) / 2;	
	nth_element(begin, mid1, end);
	return *mid + (*mid1 - *mid) / 2.0;
      }
      else { // odd number of elements
	it_t mid = begin + ( end - begin ) / 2;
	nth_element(begin, mid, end);
	return *mid;
      }
    }
  };

  // This version of median will change the order of elements.  
  template <class T> T median_unsafe(vector<T> &data) {
    //if(data.size() == 0) { std::cerr << "no data for median()" << std::endl; }
    return _median<T>()(data.begin(), data.end());
  }

  // Computes the median in O(N) time. By making a copy this version assures no re-ordering.
  template <class T> T median(vector<T> &data) {
    vector<T> data2(data.size());
    for(size_t i = 0; i < data.size(); i++) data2[i] = data[i];
    return median_unsafe<T>(data2);
  }

  template<typename T> string toString(const T& t, bool *ok = NULL) {
    ostringstream stream;
    stream << t;
    if(ok != NULL) *ok = stream.fail() == false;
    return stream.str();
  }

  template<typename T> T fromString(const string& s, bool *ok = NULL) {
    istringstream stream (s);
    T t;
    stream >> t;
    if(ok != NULL) *ok = stream.fail() == false;
    return t;
  }  

  template <typename T> T fromString(const string& s, T &defval) {
    bool ok;
    T val = fromString<T>(s, &ok);
    if(!ok) return defval; else return val;
  }

  // Priority queue adapted from Sedgewick's Algorithms in C++ Part 5.
  // Priority queue keeps track of heap positioin to vertex ID and
  // vertex ID to heap position mappings.
  template <typename T> struct PQi {
    vector<T> &a;    // vertex ID to wt or distance
    vector<int> p2v; // heap position -> vertex ID
    vector<int> v2p; // vertex ID -> heap position
    int N;
    void exch(int i, int j) { 
      int t = p2v[i]; 
      p2v[i] = p2v[j]; 
      p2v[j] = t;
      v2p[p2v[i]] = i; 
      v2p[p2v[j]] = j; 
    }
    void heapifyUp(int k) { 
      // If MIN-heap priority is violated, exchange with parent.
      while (k > 1 && a[p2v[k/2]] > a[p2v[k]]) { 
	exch(k, k/2); 
	k = k/2; // go to parent.
      } 
    }
    void heapifyDown(int k, int N) {
      while(2*k <= N) {
	int j = 2 * k; // go to child
        if (j < N && a[p2v[j]] > a[p2v[j+1]]) j++; // pick smaller child
	// Is the heap property restored, parent is not greater than
	// the smallest child.
        if (!(a[p2v[k]] > a[p2v[j]])) break;
	exch(k,j);
	k = j;
      }
    }
    PQi(int N, vector<T> &a) : a(a), p2v(N+1, 0), v2p(N+1, 0), N(0) { }
    bool empty() { return N == 0; }
    void insert(int v) { 
      p2v[++N] = v; 
      v2p[v] = N;   
      heapifyUp(N); 
    }
    // Extract the minimum and fix the heap property.
    int extractMin() { 
      exch(1, N); 
      heapifyDown(1, N-1); 
      return p2v[N--]; 
    }
    // Item k's value decreased.
    void decreaseKey(int v) { heapifyUp(v2p[v]); }
  };

};

#endif 
