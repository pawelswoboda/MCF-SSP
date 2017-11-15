#include "catch.hpp"
#include "../mcmf.hxx"

using namespace CS2_CPP;

TEST_CASE( "test problem", "[test problem]" ) {
	const int num_nodes = 6;
	const int num_arcs = 8;
	MCMF_CS2<> mcf( num_nodes, num_arcs);

   mcf.set_arc( 0, 1, 0, 4, 1);
   mcf.set_arc( 0, 2, 0, 8, 5);
   mcf.set_arc( 1, 2, 0, 5, 0);
   mcf.set_arc( 2, 4, 0, 10, 1);
   mcf.set_arc( 3, 1, 0, 8, 1);
   mcf.set_arc( 3, 5, 0, 8, 1);
   mcf.set_arc( 4, 3, 0, 8, 0);
   mcf.set_arc( 4, 5, 0, 8, 9);
   mcf.set_excess( 0, 10);
   mcf.set_excess( 5, -10);

   SECTION("solve mcf") {
      long obj = mcf.run_cs2();
      REQUIRE(obj == 70);
      // after solving, arcs are reordered lexicographically. Note that reverse arc is implicitly added by set_arc as well.
      for(int e=0; e<16; ++e) {
         std::cout << mcf.get_arc_tail(e) << " -> " << mcf.get_arc_head(e) << "; cost = " << mcf.get_cost(e) << "\n";
      }
      //REQUIRE(mcf.compute_objective_cost() == 70);

      // correct arc ordering
      REQUIRE(mcf.get_arc_tail(0) == 0); REQUIRE(mcf.get_arc_head(0) == 1);
      REQUIRE(mcf.get_arc_tail(1) == 0); REQUIRE(mcf.get_arc_head(1) == 2);
      REQUIRE(mcf.get_arc_tail(2) == 1); REQUIRE(mcf.get_arc_head(2) == 0);
      REQUIRE(mcf.get_arc_tail(3) == 1); REQUIRE(mcf.get_arc_head(3) == 2);
      REQUIRE(mcf.get_arc_tail(4) == 1); REQUIRE(mcf.get_arc_head(4) == 3);
      REQUIRE(mcf.get_arc_tail(5) == 2); REQUIRE(mcf.get_arc_head(5) == 0);
      REQUIRE(mcf.get_arc_tail(6) == 2); REQUIRE(mcf.get_arc_head(6) == 1);
      REQUIRE(mcf.get_arc_tail(7) == 2); REQUIRE(mcf.get_arc_head(7) == 4);
      REQUIRE(mcf.get_arc_tail(8) == 3); REQUIRE(mcf.get_arc_head(8) == 1);
      REQUIRE(mcf.get_arc_tail(9) == 3); REQUIRE(mcf.get_arc_head(9) == 4);
      REQUIRE(mcf.get_arc_tail(10) == 3); REQUIRE(mcf.get_arc_head(10) == 5);
      REQUIRE(mcf.get_arc_tail(11) == 4); REQUIRE(mcf.get_arc_head(11) == 2);
      REQUIRE(mcf.get_arc_tail(12) == 4); REQUIRE(mcf.get_arc_head(12) == 3);
      REQUIRE(mcf.get_arc_tail(13) == 4); REQUIRE(mcf.get_arc_head(13) == 5);
      REQUIRE(mcf.get_arc_tail(14) == 5); REQUIRE(mcf.get_arc_head(14) == 3);
      REQUIRE(mcf.get_arc_tail(15) == 5); REQUIRE(mcf.get_arc_head(15) == 4);

      REQUIRE(mcf.starting_arc(0) == 0); REQUIRE(mcf.no_arcs(0) == 2);
      REQUIRE(mcf.starting_arc(1) == 2); REQUIRE(mcf.no_arcs(1) == 3);
      REQUIRE(mcf.starting_arc(2) == 5); REQUIRE(mcf.no_arcs(2) == 3);
      REQUIRE(mcf.starting_arc(3) == 8); REQUIRE(mcf.no_arcs(3) == 3);
      REQUIRE(mcf.starting_arc(4) == 11); REQUIRE(mcf.no_arcs(4) == 3);
      REQUIRE(mcf.starting_arc(5) == 14); REQUIRE(mcf.no_arcs(5) == 2);

      // correct flow values
      // here I must know the original capacities
      //REQUIRE(mcf.get_flow(0) == 4);
      //REQUIRE(mcf.get_flow(1) == 6);
      //REQUIRE(mcf.get_flow(2) == 4);
      //REQUIRE(mcf.get_flow(3) == 10);
      //REQUIRE(mcf.get_flow(4) == 8);
      //REQUIRE(mcf.get_flow(5) == 0);
      //REQUIRE(mcf.get_flow(6) == 8);
      //REQUIRE(mcf.get_flow(7) == 2);

      REQUIRE(4-mcf.get_residual_flow(0) == 4);
      REQUIRE(8-mcf.get_residual_flow(1) == 6);
      REQUIRE(5-mcf.get_residual_flow(3) == 4);
      REQUIRE(10-mcf.get_residual_flow(7) == 10);
      REQUIRE(8-mcf.get_residual_flow(8) == 0);
      REQUIRE(8-mcf.get_residual_flow(10) == 8);
      REQUIRE(8-mcf.get_residual_flow(12) == 8);
      REQUIRE(8-mcf.get_residual_flow(13) == 2);

      REQUIRE(mcf.get_flow(0) == 4);
      REQUIRE(mcf.get_flow(1) == 6);
      REQUIRE(mcf.get_flow(3) == 4);
      REQUIRE(mcf.get_flow(7) == 10);
      REQUIRE(mcf.get_flow(8) == 0);
      REQUIRE(mcf.get_flow(10) == 8);
      REQUIRE(mcf.get_flow(12) == 8);
      REQUIRE(mcf.get_flow(13) == 2);

      // complementary slackness
      REQUIRE(mcf.get_reduced_cost(0) <= 0);
      REQUIRE(mcf.get_reduced_cost(1) == 0);
      REQUIRE(mcf.get_reduced_cost(3) == 0);
      REQUIRE(mcf.get_reduced_cost(7) <= 0);
      REQUIRE(mcf.get_reduced_cost(8) >= 0);
      REQUIRE(mcf.get_reduced_cost(10) <= 0);
      REQUIRE(mcf.get_reduced_cost(12) <= 0);
      REQUIRE(mcf.get_reduced_cost(13) == 0);
   }
}



