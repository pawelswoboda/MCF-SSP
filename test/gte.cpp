#include "../mcf_ssp.hxx"
#include "instances.h"
#include "test.h"

using namespace MCF;

int main()
{
  for(auto e : gte) {
    std::cout << "testing " << e.file << "\n";
    test_instance<long>(e.file, e.objective);
    test_instance<double>(e.file, e.objective);
  }
}
