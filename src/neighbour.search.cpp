#include <Rcpp.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <unordered_map>

using namespace std;
using namespace Rcpp;

//-----------------------------------------------
// Exact Match
//-----------------------------------------------


// if patterns is emptry, than use input sequences in trie for search.
// [[Rcpp::export(".exact_search")]]
std::vector<int> exact_search(const std::vector<std::string>& vec, const std::vector<std::string>& patterns, int max_error = 1, bool verbose = true) {
  std::vector<int> res;
  res.reserve(patterns.size() * 4);
  unordered_multimap<std::string, int> string_set;
  for (int i = 0; i < vec.size(); i++) {
    string_set.insert(std::pair<std::string, int>(vec[i], i));
  }
  for (int j = 0; j < patterns.size(); j++) {
    std::pair<unordered_map<std::string, int>::iterator,
              unordered_map<std::string, int>::iterator> i = string_set.equal_range(patterns[j]);
    if (i.first != string_set.end()) {
      for (unordered_map<std::string, int>::iterator k = i.first; k != i.second; k++) {
        res.push_back((*k).second + 1);
        res.push_back(j + 1);
      }
    }
  }
  return res;
}


// [[Rcpp::export(".exact_search_list")]]
List exact_search_list(const std::vector<std::string>& vec, const List patterns_list, int max_error = 1, bool verbose = true) {
  List out = List();
  for (int list_ind = 0; list_ind < patterns_list.length(); list_ind++) {
    std::vector<int> res;
    std::vector<std::string> patterns = as<std::vector<std::string> >(out[list_ind]);
    res.reserve(patterns.size() * 4);
    unordered_multimap<std::string, int> string_set;
    for (int i = 0; i < vec.size(); i++) {
      string_set.insert(std::pair<std::string, int>(vec[i], i));
    }
    for (int j = 0; j < patterns.size(); j++) {
      std::pair<unordered_map<std::string, int>::iterator,
                unordered_map<std::string, int>::iterator> i = string_set.equal_range(patterns[j]);
      if (i.first != string_set.end()) {
        for (unordered_map<std::string, int>::iterator k = i.first; k != i.second; k++) {
          res.push_back((*k).second + 1);
          res.push_back(j + 1);
        }
      }
    }
    out[list_ind] = 0;
  }
  return out;
}


//-----------------------------------------------
// Hamming Distance
//-----------------------------------------------


bool hamming_distance_check(const std::string& alpha, const std::string& beta, int max_error = 1) {
  if (alpha.size() != beta.size()) {
    return false;
  }

  int err = 0;
  for (int i = 0; i < alpha.size(); i++) {
    err += alpha[i] != beta[i];
    if (err > max_error) {
      return false;
    }
  }
  
  return true;
}

// if patterns is empty, than use input sequences in trie for search.
// [[Rcpp::export(".hamming_search")]]
std::vector<int> hamming_search(const std::vector<std::string>& vec, const std::vector<std::string>& patterns, int max_error = 1, bool verbose = true) {
  std::vector<int> res;
  res.reserve(patterns.size() * 4);
  for (int i = 0; i < vec.size(); i++) {
    for (int j = 0; j < patterns.size(); j++) {
      if (hamming_distance_check(vec[i], patterns[j], max_error)) {
        res.push_back(i + 1);
        res.push_back(j + 1);
      }
    }
  }
  return res;
}


//-----------------------------------------------
// Levenshtein Distance
//-----------------------------------------------


#define MAGIC_NUMBER 27

struct trie
{
    struct nucmap {

        nucmap() {
            _data = new trie*[MAGIC_NUMBER];
            for (int i = 0; i < MAGIC_NUMBER; i++) {
                _data[i] = NULL;
            }
        }

        ~nucmap() {
            for (int i = 0; i < MAGIC_NUMBER; i++) {
                delete _data[i];
            }
            delete [] _data;
        }

        trie* operator[](char letter) {
            return _data[letter - 'A'];
        }

        trie* addTrie(char letter) {
            _data[letter - 'A'] = new trie();
            return _data[letter - 'A'];
        }

        trie **_data;
    };

    typedef nucmap next_t;

    // The set with all the letters which this node is prefix
    next_t next;

    //int index;
    vector<int> index;

    trie() : next(nucmap()) { 
      //index = 0;
      index.reserve(2);
    }
    
    ~trie() {}
    

    void insert(string w, int w_index = 0)
    {
  w = '[' + w;

        trie* n = this;
        for (int i = 0; i < w.size(); ++i) {
            if (!n->next[w[i]]) {
                n->next.addTrie(w[i]);
            }
            n = n->next[w[i]];
        }

        //n->index = w_index;
        n->index.push_back(w_index);
    }
};


//
std::vector<int> search_impl(trie* tree, char ch, int *last_row, int sz, const string& word, int min_cost)
{
    int *current_row = new int[sz + 1];
    current_row[0] = last_row[0] + 1;

    // Calculate the min cost of insertion, deletion, match or substution
    int insert_or_del, replace;
    for (int i = 1; i < sz + 1; ++i) {
        insert_or_del = min(current_row[i-1] + 1, last_row[i] + 1);
        replace = (word[i-1] == ch) ? last_row[i-1] : (last_row[i-1] + 1);

        current_row[i] = min(insert_or_del, replace);
    }

    // When we find a cost that is less than the min_cost, is because
    // it is the minimum until the current row, so we update
    std::vector<int> res;
    //if ((current_row[sz] < min_cost) && (tree->index)) {
    if ((current_row[sz] < min_cost) && (tree->index.size() != 0)) {
        //res.push_back(tree->index);
        res.insert(res.end(), tree->index.begin(), tree->index.end());
    }

    // If there is an element wich is smaller than the current minimum cost,
    //  we can have another cost smaller than the current minimum cost
    if (*min_element(current_row, current_row + sz + 1) < min_cost) {
        for (int i = 'A'; i < 'A' + MAGIC_NUMBER; ++i) {
            if (tree->next[i]) {
                std::vector<int> tmp = search_impl(tree->next[i], i, current_row, sz, word, min_cost);
                if (tmp.size() > 0) {
                  res.insert(res.end(), tmp.begin(), tmp.end());
                }
            }
        }
    }

    delete [] current_row;

    return res;
}

std::vector<int> search(string word, int min_cost, trie* tree)
{
  word = '[' + word;

    int sz = word.size();

    int *current_row = new int[sz + 1];

    // Naive DP initialization
    for (int i = 0; i < sz + 1; ++i) current_row[i] = i;


    std::vector<int> res;
    // For each letter in the root map wich matches with a
    //  letter in word, we must call the search
    for (int i = 0 ; i < sz; ++i) {
        if (tree->next[word[i]]) {
            std::vector<int> tmp = search_impl(tree->next[word[i]], word[i], current_row, sz, word, min_cost);
            if (tmp.size() > 0) {
              res.insert(res.end(), tmp.begin(), tmp.end());
            }
        }
    }

    delete [] current_row;

    return res;
}

// if patterns is emptry, than use input sequences in trie for search.
// [[Rcpp::export(".levenshtein_search")]]
std::vector<int> levenshtein_search(const std::vector<std::string>& vec, const std::vector<std::string>& patterns, int max_error = 1, bool verbose = true)
{
        // The tree
    trie tree;
    
    // The minimum cost of a given word to be changed to a word of the dictionary
    int min_cost = max_error + 1;

    for (int i = 0; i < vec.size(); i++) {
        tree.insert(vec[i], i + 1);
    }

    vector<int> res;
    res.reserve(patterns.size() * 4);
    for (int i = 0; i < patterns.size(); i++) {
//      if (verbose && i % 100000 == 25) {
//        cout << i << "/" << patterns.size() << endl;
//      }
      vector<int> tmp = search(patterns[i], min_cost, &tree);
      for (int j = 0; j < tmp.size(); j++) {
        res.push_back(tmp[j]);
        res.push_back(i + 1);
      }
    }
    return res;
}