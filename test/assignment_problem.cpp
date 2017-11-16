#include "../mcf_ssp.hxx"
#include "test.h"
#include <random>
#include <vector>

using namespace MCF; 

int main()
{
   std::random_device rd;     // only used once to initialise (seed) engine
   std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
   std::uniform_int_distribution<long> uni(0,10); // guaranteed unbiased

   for(int run =0; run<100; ++run) {
      const int num_nodes = 6;
      const int num_arcs = 9;
      SSP<long,long> mcf( num_nodes, num_arcs);

      for(int i=0; i<3; ++i) {
         for(int j=0; j<3; ++j) {
            mcf.add_edge(i,3+j,0,1,uni(rng));
         }
      }
      for(int i=0; i<3; ++i) {
         mcf.add_node_excess(i,1);
         mcf.add_node_excess(3+i,-1);
      }
      mcf.order();
      mcf.solve();

      // read out solutions
      std::vector<int> flow(9);
      for(int i=0; i<9; ++i) {
         flow[i] = mcf.flow(i);
      }

      /*
      // construct reverse problems. It should have flow values < num_nodes
      SSP<long,long> mcf_reverse( num_nodes, num_arcs);
      for(int i=0; i<3; ++i) {
         for(int j=0; j<3; ++j) {
            test(flow[i*3 + j] >= 0);
            test(flow[i*3 + j] <= 1);
            if(flow[i*3 + j] == 0)
               mcf_reverse.add_edge(i,3+j,1,1000,mcf.cost(i*3+j));
            else
               mcf_reverse.add_edge(3+j,i,1,1000,mcf.cost(i*3+j));
         }
      }
      for(int i=0; i<3; ++i) {
         mcf_reverse.add_node_excess(i,0);
         mcf_reverse.add_node_excess(3+i,0);
      }
      mcf_reverse.solve();
      for(int i=0; i<num_arcs; ++i) {
         test(mcf_reverse.flow(i) < 3);
      }
      */
   }
}
