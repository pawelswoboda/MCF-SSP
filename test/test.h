#include "../mcf_ssp.hxx"
#include <stdexcept>
#include <string>

inline void test(const bool& pred)
{
  if (!pred)
    throw std::runtime_error("Test failed.");
}

template<typename CostType>
void test_instance(const std::string& filename, const long mcf_cost)
{
  auto* f = MCF::read_dimacs_file<int,CostType>(filename);
  f->order();
  auto objective = f->solve();
  std::cout << "objective value = " << objective << "\n";
  //f->print_flow();
  test(f->TestOptimality());
  test(objective == f->objective());
  test(objective == mcf_cost);
  //test(f->TestCosts());
  delete f;
}

