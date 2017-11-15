#include "catch.hpp"
#include "../mcmf.hxx"
#include <random>


using namespace CS2_CPP;

TEST_CASE( "dynamic problem", "[dynamic]" ) {
   std::random_device rd;     // only used once to initialise (seed) engine
   std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
   std::uniform_int_distribution<long> uni(0,10000); // guaranteed unbiased
   // build assignment problems and change cost after having found a solution

   /*
   SECTION( "cost update" ) {
      for(int run =0; run<100; ++run) {
         const int num_nodes = 6;
         const int num_arcs = 9;
         MCMF_CS2_STAT<> mcf( num_nodes, num_arcs);

         for(int i=0; i<3; ++i) {
            for(int j=0; j<3; ++j) {
               mcf.set_arc(i,3+j,0,1,uni(rng));
            }
         }
         for(int i=0; i<3; ++i) {
            mcf.set_excess(i,1);
            mcf.set_excess(3+i,-1);
         }
         long orig_cost = mcf.run_cs2();

         // read out solutions
         std::vector<int> flow(num_arcs);
         for(int i=0; i<num_arcs; ++i) {
            flow[i] = mcf.get_flow(i);
         }

         REQUIRE(flow[0] + flow[1] + flow[2] == 1);
         for(int i=0; i<3; ++i) {
            if(flow[i] == 1) { // increase cost by one and check whether cost of optimal solution is also increased by one
               mcf.update_cost(i, mcf.get_cost(i) + 2);
            }
         }
         long new_cost = mcf.cs2_cost_restart();
         REQUIRE(orig_cost <= new_cost);
         REQUIRE(new_cost <= orig_cost + 2);
      }
   }
   */

   SECTION( "capacity update" ) {
      for(int run =0; run<100; ++run) {
         const int num_nodes = 6;
         const int num_arcs = 9;
         MCMF_CS2_STAT<> mcf( num_nodes, num_arcs);

         for(int i=0; i<3; ++i) {
            for(int j=0; j<3; ++j) {
               mcf.set_arc(i,3+j,0,1,uni(rng));
            }
         }
         for(int i=0; i<3; ++i) {
            mcf.set_excess(i,1);
            mcf.set_excess(3+i,-1);
         }
         long orig_cost = mcf.run_cs2();
         std::vector<int> flow(num_arcs);
         for(int i=0; i<num_arcs; ++i) {
            flow[i] = mcf.get_flow(i);
         }


         for(int i=0; i<3; ++i) {
            for(int j=0; j<3; ++j) {
               mcf.set_cap(3*i + j, 2);
            }
         }
         for(int i=0; i<3; ++i) {
            mcf.inc_excess(i,1);
            mcf.inc_excess(3+i,-1);
         }

         long new_cost = mcf.cs2_restart();
         for(int i=0; i<num_arcs; ++i) {
            REQUIRE(2*flow[i] == mcf.get_flow(i));
         }
         auto nc = mcf.compute_objective_cost();
         assert(nc == new_cost);
         assert(2*orig_cost == new_cost);
         REQUIRE(2*orig_cost == new_cost);
      }
   }
}
