#include "CME212/Util.hpp"

#include "Graph.hpp"


static unsigned fail_count = 0;

template <typename T>
void sf_print(T a, std::string msg = "") {
  (void) a;
  std::cerr << msg << " [Success]" << std::endl;
}

void sf_print(bool sf, std::string msg = "") {
  if (sf)
    std::cerr << msg << " [Success]" << std::endl;
  else {
    std::cerr << msg << " [FAIL]" << std::endl;
    ++fail_count;
  }
}


int main()
{
  using GraphType = Graph<int, int>;
  using Node = GraphType::node_type;

  GraphType g;

  Point p(CME212::random(), CME212::random(), CME212::random());

  std::cerr << "Node is " << sizeof(Node) << " bytes";
  sf_print(sizeof(Node) <= 16);

  sf_print(g.add_node(p), "Inserting Node");

  sf_print(g.num_nodes() == 1, "Graph has 1 Node");

  sf_print(g.node(0), "Getting Node");
  Node node = g.node(0);

  sf_print(node.position() == p, "Node position conserved");
  sf_print(node.index() == 0, "Index is 0");
  sf_print(node.value() == 0, "Node value is default");
  std::cerr << "Setting node value...";
  g.node(0).value() = 2;
  sf_print(node.value() == 2);

  g.remove_node(node);
  sf_print(g.num_nodes() == 0, "Removing Node...");

  sf_print(g.num_nodes() == 0, "Graph has 0 Nodes");

  std::cerr << "Adding 100 Nodes...";
  for (int k = 0; k < 100; ++k)
    g.add_node(Point(CME212::random(), CME212::random(), CME212::random()), k);
  sf_print(g.num_nodes() == 100);
  sf_print(g.node(0).value() == 0, "Node value set (0)");
  sf_print(g.node(99).value() == 99, "Node value set (99)");

  // Remove 50 Nodes...
  for (unsigned k = 0; k < 50; ++k) {
    unsigned n = unsigned(CME212::random(0, g.num_nodes()));
    g.remove_node(g.node(n));
  }

  sf_print(g.num_nodes() == 50, "Removed 50 nodes");

  std::cerr << "Checking Node indices...";
  bool succ = true;
  for (unsigned k = 0; succ && k < 50; ++k) {
    Node node = g.node(k);
    if (node.index() != k)
      succ = false;
  }
  sf_print(succ);

  std::cerr << "Clearing...";
  g.clear();
  sf_print(g.num_nodes() == 0, "Graph has 0 Nodes");

  GraphType g2;
  std::cerr << "Adding 50 Nodes to Graph1 and Graph2...";
  for (unsigned k = 0; k < 50; ++k) {
    Point p(CME212::random(), CME212::random(), CME212::random());
    g.add_node(p);
    g2.add_node(p);
  }
  sf_print(g.num_nodes() == 50 && g2.num_nodes() == 50);

  sf_print(g2.node(23) == g2.node(23), "G2-G2 Node comparison ==");
  sf_print(g2.node(23) != g2.node(21), "G2-G2 Node comparison !=");
  sf_print(g2.node(23) != g.node(23), "G2-G1 Node comparison !=");
  sf_print(g2.node(23) != g.node(21), "G2-G1 Node comparison !=");

  sf_print(g.node(23) < g.node(21) || g.node(23) > g.node(21),
           "G1-G1 Node comparison < >");
  sf_print(g.node(23) < g2.node(21) || g.node(23) > g2.node(21),
           "G1-G2 Node comparison < >");

  if (fail_count) {
    std::cerr << "\n" << fail_count
	      << (fail_count > 1 ? " FAILURES" : " FAILURE") << std::endl;
    return 1;
  } else
    return 0;
}
