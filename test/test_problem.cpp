#include "test.h"
#include "../mcf_ssp.hxx"

using namespace MCF;

int main()
{
	const int num_nodes = 6;
	const int num_arcs = 8;
	SSP<long,long> mcf( num_nodes, num_arcs);

   mcf.add_edge( 0, 1, 0, 4, 1);
   mcf.add_edge( 0, 2, 0, 8, 5);
   mcf.add_edge( 1, 2, 0, 5, 0);
   mcf.add_edge( 2, 4, 0, 10, 1);
   mcf.add_edge( 3, 1, 0, 8, 1);
   mcf.add_edge( 3, 5, 0, 8, 1);
   mcf.add_edge( 4, 3, 0, 8, 0);
   mcf.add_edge( 4, 5, 0, 8, 9);
   mcf.add_node_excess( 0, 10);
   mcf.add_node_excess( 5, -10);

   mcf.order();

   // test copy constructor
   SSP<long,long> mcf2(mcf);

   auto obj = mcf.solve();
   test(obj == 70);
   // after solving, arcs are reordered lexicographically. Note that reverse arc is implicitly added by add_edge as well.
   for(int e=0; e<16; ++e) {
     std::cout << mcf.tail(e) << " -> " << mcf.head(e) << "; cost = " << mcf.cost(e) << "\n";
   }

   // correct arc ordering
   test(mcf.tail(0) == 0); test(mcf.head(0) == 1); test(mcf.lower_bound(0) == 0); test(mcf.upper_bound(0) == 4); 
   test(mcf.tail(1) == 0); test(mcf.head(1) == 2); test(mcf.lower_bound(1) == 0); test(mcf.upper_bound(1) == 8); 
   test(mcf.tail(2) == 1); test(mcf.head(2) == 0); test(mcf.lower_bound(2) == 4); test(mcf.upper_bound(2) == 0); 
   test(mcf.tail(3) == 1); test(mcf.head(3) == 2); test(mcf.lower_bound(3) == 0); test(mcf.upper_bound(3) == 5); 
   test(mcf.tail(4) == 1); test(mcf.head(4) == 3); test(mcf.lower_bound(4) == 8); test(mcf.upper_bound(4) == 0); 
   test(mcf.tail(5) == 2); test(mcf.head(5) == 0); test(mcf.lower_bound(5) == 8); test(mcf.upper_bound(5) == 0); 
   test(mcf.tail(6) == 2); test(mcf.head(6) == 1); test(mcf.lower_bound(6) == 5); test(mcf.upper_bound(6) == 0); 
   test(mcf.tail(7) == 2); test(mcf.head(7) == 4); test(mcf.lower_bound(7) == 0); test(mcf.upper_bound(7) == 10); 
   test(mcf.tail(8) == 3); test(mcf.head(8) == 1); test(mcf.lower_bound(8) == 0); test(mcf.upper_bound(8) == 8); 
   test(mcf.tail(9) == 3); test(mcf.head(9) == 4); test(mcf.lower_bound(9) == 8); test(mcf.upper_bound(9) == 0); 
   test(mcf.tail(10) == 3); test(mcf.head(10) == 5); test(mcf.lower_bound(10) == 0); test(mcf.upper_bound(10) == 8); 
   test(mcf.tail(11) == 4); test(mcf.head(11) == 2); test(mcf.lower_bound(11) == 10); test(mcf.upper_bound(11) == 0); 
   test(mcf.tail(12) == 4); test(mcf.head(12) == 3); test(mcf.lower_bound(12) == 0); test(mcf.upper_bound(12) == 8); 
   test(mcf.tail(13) == 4); test(mcf.head(13) == 5); test(mcf.lower_bound(13) == 0); test(mcf.upper_bound(13) == 8); 
   test(mcf.tail(14) == 5); test(mcf.head(14) == 3); test(mcf.lower_bound(14) == 8); test(mcf.upper_bound(14) == 0); 
   test(mcf.tail(15) == 5); test(mcf.head(15) == 4); test(mcf.lower_bound(15) == 8); test(mcf.upper_bound(15) == 0); 

   test(mcf.first_outgoing_arc(0) == 0); test(mcf.no_outgoing_arcs(0) == 2);
   test(mcf.first_outgoing_arc(1) == 2); test(mcf.no_outgoing_arcs(1) == 3);
   test(mcf.first_outgoing_arc(2) == 5); test(mcf.no_outgoing_arcs(2) == 3);
   test(mcf.first_outgoing_arc(3) == 8); test(mcf.no_outgoing_arcs(3) == 3);
   test(mcf.first_outgoing_arc(4) == 11); test(mcf.no_outgoing_arcs(4) == 3);
   test(mcf.first_outgoing_arc(5) == 14); test(mcf.no_outgoing_arcs(5) == 2);

   // correct flow values
   // here I must know the original capacities
   //test(mcf.flow(0) == 4);
   //test(mcf.flow(1) == 6);
   //test(mcf.flow(2) == 4);
   //test(mcf.flow(3) == 10);
   //test(mcf.flow(4) == 8);
   //test(mcf.flow(5) == 0);
   //test(mcf.flow(6) == 8);
   //test(mcf.flow(7) == 2);

   test(4-mcf.residual_capacity(0) == 4);
   test(8-mcf.residual_capacity(1) == 6);
   test(5-mcf.residual_capacity(3) == 4);
   test(10-mcf.residual_capacity(7) == 10);
   test(8-mcf.residual_capacity(8) == 0);
   test(8-mcf.residual_capacity(10) == 8);
   test(8-mcf.residual_capacity(12) == 8);
   test(8-mcf.residual_capacity(13) == 2);

   test(mcf.flow(0) == 4);
   test(mcf.flow(1) == 6);
   test(mcf.flow(3) == 4);
   test(mcf.flow(7) == 10);
   test(mcf.flow(8) == 0);
   test(mcf.flow(10) == 8);
   test(mcf.flow(12) == 8);
   test(mcf.flow(13) == 2);

   // complementary slackness
   test(mcf.reduced_cost(0) <= 0);
   test(mcf.reduced_cost(1) == 0);
   test(mcf.reduced_cost(3) == 0);
   test(mcf.reduced_cost(7) <= 0);
   test(mcf.reduced_cost(8) >= 0);
   test(mcf.reduced_cost(10) <= 0);
   test(mcf.reduced_cost(12) <= 0);
   test(mcf.reduced_cost(13) == 0);

   mcf2.order();
   mcf2.print_flow();
   mcf2.solve();
   test(mcf2.objective() == obj);

   SSP<long,long> mcf3;
   mcf3 = mcf;
   test(mcf3.objective() == obj);

   SSP<long,long> mcf4;
   mcf4 = std::move(mcf);
   test(mcf4.objective() == obj);
   test(mcf.objective() == 0);
}
