#include "../ssp_mcf.hxx"
#include "test_instances.h"
#include "test.h"

using namespace MCF;

int main()
{
  // SSP is not able to solve those in a reasonable time
  for(auto e : netgen) {
    std::cout << "testing " << e.file << "\n";
    test_instance(e.file, e.objective);
  }
}


