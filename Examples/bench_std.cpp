#include <iostream>
#include <algorithm>

#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>

#include "CME212/Util.hpp"

/** Quickie RAII Timer */
struct Timer {
  std::string msg;
  CME212::Clock clock;

  Timer(const std::string& _msg)
      : msg(_msg) {
    clock.start();
  }

  ~Timer() {
    double elapsed = clock.seconds();
    std::cout << msg << ": " << elapsed << "s" << std::endl;
  }
};


int main()
{
  std::vector<int> a;
  std::set<int> s;
  std::map<int,int> m;
  std::unordered_set<int> us;
  std::unordered_map<int,int> um;

  int N = 10;
  std::cout << "N = " << N << std::endl;

  { Timer timer("Vector push_back");
    for (int k = 0; k < N; ++k)
      a.push_back(k);
  }

  { Timer timer("Set insert");
    for (int k = 0; k < N; ++k)
      s.insert(k);
  }

  { Timer timer("Map insert");
    for (int k = 0; k < N; ++k)
      m[k] = k;
  }

  { Timer timer("Unordered set insert");
    for (int k = 0; k < N; ++k)
      us.insert(k);
  }

  { Timer timer("Unordered map insert");
    for (int k = 0; k < N; ++k)
      um[k] = k;
  }

  int R = 100;

  { Timer timer("Vector iterator traversal");
    for (int r = 0; r < R; ++r)
      for (auto it = a.begin(); it != a.end(); ++it)
        assert(*it >= 0);
  }

  { Timer timer("Set traversal");
    for (int r = 0; r < R; ++r)
      for (auto it = s.begin(); it != s.end(); ++it)
        assert(*it >= 0);
  }

  { Timer timer("Map traversal");
    for (int r = 0; r < R; ++r)
      for (auto it = m.begin(); it != m.end(); ++it)
        assert(it->second >= 0);
  }

  { Timer timer("Unordered set traversal");
    for (int r = 0; r < R; ++r)
      for (auto it = us.begin(); it != us.end(); ++it)
        assert(*it >= 0);
  }

  { Timer timer("Unordered map traversal");
    for (int r = 0; r < R; ++r)
      for (auto it = um.begin(); it != um.end(); ++it)
        assert(it->second >= 0);
  }

  int S = 10000;

  { Timer timer("Vector op[]");
    for (int k = 0; k < S; ++k)
      assert(a[int(CME212::random(0,N))] >= 0);
  }

  { Timer timer("Set.find");
    for (int k = 0; k < S; ++k)
      assert(s.find(int(CME212::random(0,N))) != s.end());
  }

  { Timer timer("Map op[]");
    for (int k = 0; k < S; ++k)
      assert(m[int(CME212::random(0,N))] >= 0);
  }

  { Timer timer("Unordered set find");
    for (int k = 0; k < S; ++k)
      assert(us.find(int(CME212::random(0,N))) != us.end());
  }

  { Timer timer("Unordered map op[]");
    for (int k = 0; k < S; ++k)
      assert(um[int(CME212::random(0,N))] >= 0);
  }

  { Timer timer("Vector std::binary_search");
    for (int k = 0; k < S; ++k)
      assert(std::binary_search(a.begin(), a.end(), int(CME212::random(0,N))));
  }

  { Timer timer("Vector std::find");
    for (int k = 0; k < S; ++k)
      assert(std::find(a.begin(), a.end(), int(CME212::random(0,N))) != a.end());
  }

  return 0;
}
