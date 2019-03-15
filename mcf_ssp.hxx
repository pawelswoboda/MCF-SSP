#ifndef MCF_SSP_HXX
#define MCF_SSP_HXX

#include <cmath>
#include <cstring>
#include <assert.h>

#include <algorithm>
#include <numeric>
#include <memory>
#include <string>

#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>


namespace MCF {

    template <typename FlowType, typename CostType> class SSP
    {
        public:
            typedef std::size_t NodeId;
            typedef std::size_t EdgeId;

            SSP();
            SSP(std::size_t NodeNum, std::size_t edgeNumMax);
            SSP(const SSP& o);
            SSP(SSP&& o);
            SSP& operator=(SSP& o);
            SSP& operator=(SSP&& o);

            template<typename _FlowType, typename _CostType>
                friend void swap(SSP<_FlowType,_CostType>&,SSP<_FlowType,_CostType>&);

            void copy_node(const SSP& o, NodeId i);
            void copy_arc(const SSP& o, EdgeId i);

            // Destructor
            ~SSP();

            void add_node_excess(NodeId i, FlowType excess);

            // first call returns 0, second 1, and so on.
            // cap, rev_cap must be non-negative. 
            // cost can be negative.
            // EdgeIds only stay unchanged when arcs are not reordered
            EdgeId add_edge(NodeId i, NodeId j, FlowType lower, FlowType upper, CostType cost);

            CostType solve();
            CostType objective() const;

            ///////////////////////////////////////////////////

            FlowType GetRCap(EdgeId e);
            void SetRCap(EdgeId e, FlowType new_rcap);
            FlowType GetReverseRCap(EdgeId e);
            void SetReverseRCap(EdgeId e, FlowType new_rcap);
            void PushFlow(EdgeId e, FlowType delta);
            void update_cost(EdgeId e, CostType delta);
            void reset_costs();

            // query functions 
            void order() { order_inter_nodes(); order_intra_nodes(); } // reorder arcs so that outgoing ones from given node are ordered consecutively
            NodeId no_nodes() const;
            EdgeId no_edges() const;
            EdgeId no_arcs() const;
            FlowType flow(const NodeId i, const EdgeId e) const; // get the flow of the e-th edge outgoing out of i
            FlowType flow(const EdgeId e) const; // get the flow of the e-th edge outgoing out of i
            CostType cost(const EdgeId e) const { assert(e >= 0 && e < 2*edgeNum); return arcs[e].cost; }
            CostType reduced_cost(const EdgeId e) const { assert(e >= 0 && e < 2*edgeNum); return arcs[e].GetRCost(); }
            CostType residual_capacity(const EdgeId e) const { assert(e >= 0 && e < 2*edgeNum); return arcs[e].r_cap; }
            NodeId tail(EdgeId e) const { return arcs[e].sister->head - nodes; }
            NodeId head(EdgeId e) const { return arcs[e].head - nodes; }
            EdgeId first_outgoing_arc(NodeId i) const;
            std::size_t no_outgoing_arcs(NodeId i) const;
            FlowType upper_bound(EdgeId i) const { assert(arc_valid(&arcs[i])); return capacity[i]; }
            FlowType lower_bound(EdgeId i) const { assert(arc_valid(&arcs[i])); const EdgeId s = arcs[i].sister - arcs; assert(arc_valid(&arcs[s])); return capacity[s]; }
            CostType potential(NodeId i) const { assert(i<no_nodes()); return nodes[i].pi; }

            // debug functions
            bool TestOptimality() const; 
            bool TestCosts() const;
            void print_flow() const;

            /////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////

        private:
            // internal variables and functions

            struct Node;
            struct Arc;

            struct Node
            {
                Arc			*firstNonsaturated;
                Arc			*firstSaturated;

                Arc			*parent;
                Node		*next; // list of nodes with positive excesses

                FlowType	excess;
                CostType	pi;
                std::size_t			flag;
                union
                {
                    std::size_t		heap_ptr;
                    Node*	next_permanent;
                };
            };

            struct Arc
            {
                Node		*head;
                Arc			*prev;    // previous arc in saturated or non-saturated list
                Arc			*next;    // next arc in saturated or non-saturated list
                Arc			*sister;	// reverse arc

                FlowType	r_cap;		// residual capacity
                CostType	cost;
                CostType GetRCost() const { return cost + head->pi - sister->head->pi; }
            };

            std::size_t		nodeNum, edgeNum, edgeNumMax;
            Node	*nodes = nullptr;
            Arc		*arcs = nullptr;
            Node*	firstActive = nullptr;
            std::size_t		counter;
            CostType mcf_cost;

            FlowType* capacity = nullptr; // original capacities from which one can compute the flows


            /////////////////////////////////////////////////////////////////////////

            struct PriorityQueue
            {
                PriorityQueue();
                ~PriorityQueue();
                void Reset();
                CostType GetKey(Node* i);
                void Add(Node* i, CostType key);
                void DecreaseKey(Node* i, CostType key);
                Node* RemoveMin(CostType& key);

                private:
                struct Item
                {
                    Node*		i;
                    CostType	key;
                }* array;
                std::size_t N, arraySize;
                void Swap(std::size_t k1, std::size_t k2);
            };

            PriorityQueue queue;

            /////////////////////////////////////////////////////////////////////////

            void SetRCap(Arc* a, FlowType new_rcap);
            void PushFlow(Arc* a, FlowType delta);

            void Init();
            void DecreaseRCap(Arc* a, FlowType delta);
            void IncreaseRCap(Arc* a, FlowType delta);
            FlowType Augment(Node* start, Node* end);
            void Dijkstra(Node* start);

            bool node_valid(NodeId i) const;
            bool arc_valid(Arc* a) const;
            void exchange(Arc* const a, Arc* const b);

            void order_inter_nodes(); // reorder arcs so that the outgoing ones from each node are ordered consecutively.
            void order_intra_nodes(); // reorder arcs outgoing from each node by head node id
    };



    ///////////////////////////////////////
    // Implementation - inline functions //
    ///////////////////////////////////////

    template <typename FlowType, typename CostType> 
        inline FlowType SSP<FlowType, CostType>::flow(EdgeId _e) const
        {
            return capacity[_e] - arcs[_e].r_cap;
        }

    template <typename FlowType, typename CostType> 
        inline FlowType SSP<FlowType, CostType>::flow(NodeId _i, EdgeId _e) const
        {
            assert(false);
            EdgeId e = (nodes[_i].first() + _e) - arcs;
            return flow(e);
        }


    template <typename FlowType, typename CostType> 
        inline typename SSP<FlowType, CostType>::NodeId SSP<FlowType, CostType>::no_nodes() const
        {
            return nodeNum;
        }
    template <typename FlowType, typename CostType> 
        inline typename SSP<FlowType, CostType>::EdgeId SSP<FlowType, CostType>::no_edges() const
        {
            return edgeNum;
        }
    template <typename FlowType, typename CostType> 
        inline typename SSP<FlowType, CostType>::EdgeId SSP<FlowType, CostType>::no_arcs() const
        {
            return 2*edgeNum;
        }

    template <typename FlowType, typename CostType> 
        inline std::size_t SSP<FlowType, CostType>::no_outgoing_arcs(NodeId i) const
        {
            assert(node_valid(i));
            std::size_t n = 0;
            for (Arc* a=nodes[i].firstSaturated; a; a=a->next) { ++n; }
            for (Arc* a=nodes[i].firstNonsaturated; a; a=a->next) { ++n; }
            return n;
        }

    // only makes sense if arcs have been ordered
    template <typename FlowType, typename CostType> 
        inline typename SSP<FlowType, CostType>::EdgeId SSP<FlowType, CostType>::first_outgoing_arc(NodeId i) const
        {
            assert(node_valid(i));
            EdgeId e = std::numeric_limits<EdgeId>::max();
            for (Arc* a=nodes[i].firstSaturated; a; a=a->next) { 
                e = std::min(e, EdgeId(a-arcs));
            }
            for (Arc* a=nodes[i].firstNonsaturated; a; a=a->next) { 
                e = std::min(e, EdgeId(a-arcs));
            }
            return e;
        }

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::add_node_excess(NodeId _i, FlowType excess)
        {
            assert(_i>=0 && _i<nodeNum);
            nodes[_i].excess += excess;
            if (nodes[_i].excess > 0 && !nodes[_i].next)
            {
                nodes[_i].next = firstActive;
                firstActive = &nodes[_i];
            }
        }

    template <typename FlowType, typename CostType> 
        inline typename SSP<FlowType, CostType>::EdgeId SSP<FlowType, CostType>::add_edge(NodeId _i, NodeId _j, FlowType lower, FlowType upper, CostType cost)
        {
            assert(_i>=0 && _i<nodeNum);
            assert(_j>=0 && _j<nodeNum);
            assert(_i!=_j && edgeNum<edgeNumMax);
            assert(upper >= 0);
            assert(lower <= 0);
            assert(lower < upper);

            Arc *a = &arcs[2*edgeNum];
            Arc *a_rev = a+1;

            capacity[2*edgeNum] = upper;
            capacity[2*edgeNum+1] = lower;

            edgeNum ++;

            Node* i = nodes + _i;
            Node* j = nodes + _j;

            a -> sister = a_rev;
            a_rev -> sister = a;
            if (upper > 0)
            {
                if (i->firstNonsaturated) i->firstNonsaturated->prev = a;
                a -> next = i -> firstNonsaturated;
                i -> firstNonsaturated = a;
            }
            else
            {
                if (i->firstSaturated) i->firstSaturated->prev = a;
                a -> next = i -> firstSaturated;
                i -> firstSaturated = a;
            }
            a->prev = nullptr;
            if (lower < 0)
            {
                if (j->firstNonsaturated) j->firstNonsaturated->prev = a_rev;
                a_rev -> next = j -> firstNonsaturated;
                j -> firstNonsaturated = a_rev;
            }
            else
            {
                if (j->firstSaturated) j->firstSaturated->prev = a_rev;
                a_rev -> next = j -> firstSaturated;
                j -> firstSaturated = a_rev;
            }
            a_rev->prev = nullptr;

            a -> head = j;
            a_rev -> head = i;
            a -> r_cap = upper;
            a_rev -> r_cap = -lower;
            a -> cost = cost;
            a_rev -> cost = -cost;

            assert(arc_valid(a) && arc_valid(a_rev));

            if (a->r_cap > 0 && a->GetRCost() < 0) PushFlow(a, a->r_cap);
            if (a_rev->r_cap > 0 && a_rev->GetRCost() < 0) PushFlow(a_rev, a_rev->r_cap);

            assert(arc_valid(a) && arc_valid(a_rev));
            return edgeNum-1;
        }

    ///////////////////////////////////////
    ///////////////////////////////////////
    ///////////////////////////////////////

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::DecreaseRCap(Arc* a, FlowType delta)
        {
            a->r_cap -= delta;
            if (a->r_cap == 0)
            {
                Node* i = a->sister->head;
                if (a->next) a->next->prev = a->prev;
                if (a->prev) a->prev->next = a->next;
                else         i->firstNonsaturated = a->next;
                a->next = i->firstSaturated;
                if (a->next) a->next->prev = a;
                a->prev = nullptr;
                i->firstSaturated = a;
            }
        }

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::IncreaseRCap(Arc* a, FlowType delta)
        {
            if (a->r_cap == 0)
            {
                Node* i = a->sister->head;
                if (a->next) a->next->prev = a->prev;
                if (a->prev) a->prev->next = a->next;
                else         i->firstSaturated = a->next;
                a->next = i->firstNonsaturated;
                if (a->next) a->next->prev = a;
                a->prev = nullptr;
                i->firstNonsaturated = a;
            }
            a->r_cap += delta;
        }

    template <typename FlowType, typename CostType> 
        inline FlowType SSP<FlowType, CostType>::GetRCap(EdgeId e)
        {
            Arc* a = &arcs[2*e];
            return a->r_cap;
        }

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::SetRCap(Arc* a, FlowType new_rcap)
        {
            assert(new_rcap >= 0);
#ifdef SSP_DEBUG
            a->cap_orig += new_rcap - a->r_cap;
#endif
            if (a->r_cap == 0)
            {
                Node* i = a->sister->head;
                if (a->next) a->next->prev = a->prev;
                if (a->prev) a->prev->next = a->next;
                else         i->firstSaturated = a->next;
                a->next = i->firstNonsaturated;
                if (a->next) a->next->prev = a;
                a->prev = nullptr;
                i->firstNonsaturated = a;
            }
            a->r_cap = new_rcap;
            if (a->r_cap == 0)
            {
                Node* i = a->sister->head;
                if (a->next) a->next->prev = a->prev;
                if (a->prev) a->prev->next = a->next;
                else         i->firstNonsaturated = a->next;
                a->next = i->firstSaturated;
                if (a->next) a->next->prev = a;
                a->prev = nullptr;
                i->firstSaturated = a;
            }
        }

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::SetRCap(EdgeId e, FlowType new_rcap)
        {
            SetRCap(&arcs[2*e], new_rcap);
        }

    template <typename FlowType, typename CostType> 
        inline FlowType SSP<FlowType, CostType>::GetReverseRCap(EdgeId e)
        {
            Arc* a = &arcs[2*e+1];
            return a->r_cap;
        }

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::SetReverseRCap(EdgeId e, FlowType new_rcap)
        {
            SetRCap(&arcs[2*e+1], new_rcap);
        }

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::PushFlow(Arc* a, FlowType delta)
        {
            if (delta < 0) { a = a->sister; delta = -delta; }
            DecreaseRCap(a, delta);
            IncreaseRCap(a->sister, delta);
            a->head->excess += delta;
            a->sister->head->excess -= delta;
            mcf_cost += delta*a->cost;
            if (a->head->excess > 0 && !a->head->next)
            {
                a->head->next = firstActive;
                firstActive = a->head;
            }
        }

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::PushFlow(EdgeId e, FlowType delta)
        {
            PushFlow(&arcs[2*e], delta);
        }

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::update_cost(EdgeId e, CostType delta)
        {
            Arc* a = &arcs[e];
            mcf_cost += delta*(capacity[e]-a->r_cap);
            a->cost += delta;
            a->sister->cost = -a->cost;

            if (a->GetRCost() > 0) a = a->sister;
            if (a->r_cap > 0 && a->GetRCost() < 0) PushFlow(a, a->r_cap);
        }

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::reset_costs()
        {
           for(EdgeId e=0; e<no_arcs(); ++e) {
              update_cost(e, -cost(e));
              assert(arcs[e].cost == 0.0);
           }
           mcf_cost = 0.0;
           assert(objective() == 0.0);
        }

    ///////////////////////////////////////
    ///////////////////////////////////////
    ///////////////////////////////////////

    template <typename FlowType, typename CostType> 
        inline SSP<FlowType, CostType>::PriorityQueue::PriorityQueue()
        {
            N = 0;
            arraySize = 16;
            array = (Item*) malloc(arraySize*sizeof(Item));
        }

    template <typename FlowType, typename CostType> 
        inline SSP<FlowType, CostType>::PriorityQueue::~PriorityQueue()
        {
            if(array)
                free(array);
        }

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::PriorityQueue::Reset()
        {
            N = 0;
        }

    template <typename FlowType, typename CostType> 
        inline CostType SSP<FlowType, CostType>::PriorityQueue::GetKey(Node* i)
        {
            return array[i->heap_ptr].key;
        }

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::PriorityQueue::Swap(std::size_t k1, std::size_t k2)
        {
            Item* a = array+k1;
            Item* b = array+k2;
            a->i->heap_ptr = k2;
            b->i->heap_ptr = k1;
            Node* i = a->i;   a->i   = b->i;   b->i   = i;
            CostType key = a->key; a->key = b->key; b->key = key;
        }

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::PriorityQueue::Add(Node* i, CostType key)
        {
            if (N == arraySize)
            {
                arraySize *= 2;
                array = (Item*) realloc(array, arraySize*sizeof(Item));
            }
            std::size_t k = i->heap_ptr = N ++;
            array[k].i = i;
            array[k].key = key;
            while (k > 0)
            {
                std::size_t k2 = (k-1)/2;
                if (array[k2].key <= array[k].key) break;
                Swap(k, k2);
                k = k2;
            }
        }

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::PriorityQueue::DecreaseKey(Node* i, CostType key)
        {
            std::size_t k = i->heap_ptr;
            array[k].key = key;
            while (k > 0)
            {
                std::size_t k2 = (k-1)/2;
                if (array[k2].key <= array[k].key) break;
                Swap(k, k2);
                k = k2;
            }
        }

    template <typename FlowType, typename CostType> 
        inline typename SSP<FlowType, CostType>::Node* SSP<FlowType, CostType>::PriorityQueue::RemoveMin(CostType& key)
        {
            if (N == 0) return nullptr;

            Swap(0, N-1);
            N --;

            std::size_t k = 0;
            while ( 1 )
            {
                std::size_t k1 = 2*k + 1, k2 = k1 + 1;
                if (k1 >= N) break;
                std::size_t k_min = (k2 >= N || array[k1].key <= array[k2].key) ? k1 : k2;
                if (array[k].key <= array[k_min].key) break;
                Swap(k, k_min);
                k = k_min;
            }

            key = array[N].key;
            return array[N].i;
        }

    template <typename FlowType, typename CostType> 
        inline SSP<FlowType, CostType>::SSP()
        : nodeNum(0),
        edgeNum(0),
        edgeNumMax(0),
        counter(0),
        mcf_cost(0),
        nodes(nullptr),
        arcs(nullptr),
        capacity(nullptr),
        firstActive(nullptr)
    {}
    template <typename FlowType, typename CostType> 
        inline SSP<FlowType, CostType>::SSP(std::size_t _nodeNum, std::size_t _edgeNumMax)
        : nodeNum(_nodeNum),
        edgeNum(0),
        edgeNumMax(_edgeNumMax),
        counter(0),
        mcf_cost(0)
    {
        if(nodeNum > 0) { nodes = (Node*) malloc(nodeNum*sizeof(Node)); }
        if(edgeNumMax > 0) { arcs = (Arc*) malloc(2*edgeNumMax*sizeof(Arc)); }
        if(edgeNumMax > 0) { capacity = (FlowType*) malloc(2*edgeNumMax*sizeof(FlowType)); }
        if ((nodeNum > 0 && !nodes) || 
                (edgeNumMax > 0 && !arcs) || 
                (edgeNumMax > 0 && !capacity) )
        { throw std::bad_alloc(); }

        std::memset(nodes, 0, nodeNum*sizeof(Node));
        std::memset(arcs, 0, 2*edgeNumMax*sizeof(Arc));
        firstActive = &nodes[nodeNum];
    }

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::copy_node(const SSP& o, NodeId i)
        {
            if(o.nodes[i].firstNonsaturated != nullptr) { nodes[i].firstNonsaturated = arcs + (o.nodes[i].firstNonsaturated - o.arcs); }
            else { nodes[i].firstNonsaturated = nullptr; }
            if(o.nodes[i].firstSaturated != nullptr) { nodes[i].firstSaturated = arcs + (o.nodes[i].firstSaturated - o.arcs); }
            else { nodes[i].firstSaturated = nullptr; }
            if(o.nodes[i].parent != nullptr) { nodes[i].parent = arcs + (o.nodes[i].parent - o.arcs); }
            else { nodes[i].parent = nullptr; }
            if(o.nodes[i].next != nullptr) { nodes[i].next = nodes + (o.nodes[i].next - o.nodes); }
            else { nodes[i].next = nullptr; }

            nodes[i].excess = o.nodes[i].excess;
            nodes[i].pi = o.nodes[i].pi;
            nodes[i].flag = o.nodes[i].flag;
        }

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::copy_arc(const SSP& o, EdgeId i)
        {
            arcs[i].head = nodes + (o.arcs[i].head - o.nodes);
            arcs[i].sister = arcs + (o.arcs[i].sister - o.arcs);

            if(o.arcs[i].next != nullptr) { arcs[i].next = arcs + (o.arcs[i].next - o.arcs); }
            else { arcs[i].next = nullptr; }
            if(o.arcs[i].prev != nullptr) { arcs[i].prev = arcs + (o.arcs[i].prev - o.arcs); }
            else { arcs[i].prev = nullptr; }

            arcs[i].r_cap = o.arcs[i].r_cap;
            arcs[i].cost = o.arcs[i].cost;

            capacity[i] = o.capacity[i];
        }

    template <typename FlowType, typename CostType> 
        inline SSP<FlowType, CostType>::SSP(const SSP& o)
        : nodeNum(o.nodeNum),
        edgeNum(o.edgeNum),
        edgeNumMax(o.edgeNumMax),
        counter(o.counter),
        mcf_cost(o.mcf_cost)
    {
        nodes = (Node*) malloc(nodeNum*sizeof(Node));
        arcs = (Arc*) malloc(2*edgeNumMax*sizeof(Arc));
        capacity = (FlowType*) malloc(2*edgeNumMax*sizeof(FlowType));
        if (!nodes || !arcs || !capacity) { throw std::bad_alloc(); }

        for(NodeId i=0; i<nodeNum; ++i) { copy_node(o,i); }
        for(EdgeId i=0; i<2*edgeNum; ++i) { copy_arc(o,i); }

        firstActive = nodes + (o.firstActive - o.nodes);
    }

    template <typename FlowType, typename CostType> 
        inline void swap(SSP<FlowType,CostType>& first, SSP<FlowType,CostType>& second)
        {
            using std::swap;

            std::swap(first.nodeNum, second.nodeNum);
            std::swap(first.edgeNum, second.edgeNum);
            std::swap(first.edgeNumMax, second.edgeNumMax);
            std::swap(first.counter, second.counter);
            std::swap(first.mcf_cost, second.mcf_cost);
            std::swap(first.firstActive, second.firstActive);

            std::swap(first.nodes, second.nodes);
            std::swap(first.arcs, second.arcs);
            std::swap(first.capacity, second.capacity);
        }

    template <typename FlowType, typename CostType> 
        inline SSP<FlowType, CostType>::SSP(SSP&& o)
        {
            swap(*this, o);
        } 

    template <typename FlowType, typename CostType> 
        inline SSP<FlowType,CostType>& SSP<FlowType, CostType>::operator=(SSP<FlowType,CostType>& o)
        {
            SSP<FlowType,CostType> o2(o);
            swap(*this, o2);
            return *this;
        }

    template <typename FlowType, typename CostType> 
        inline SSP<FlowType,CostType>& SSP<FlowType, CostType>::operator=(SSP<FlowType,CostType>&& o)
        {
            swap(*this, o);
            return *this; 
        }

    template <typename FlowType, typename CostType> 
        inline SSP<FlowType, CostType>::~SSP()
        {
            if(nodes != nullptr) free(nodes);
            if(arcs != nullptr) free(arcs);
            if(capacity != nullptr) free(capacity);
        }

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::Init()
        {
            Node* i;
            Arc* a;

            for (a=arcs; a<arcs+2*edgeNum; a++)
            {
                if (a->r_cap > 0 && a->GetRCost() < 0) PushFlow(a, a->r_cap);
            }

            Node** lastActivePtr = &firstActive;
            for (i=nodes; i<nodes+nodeNum; i++)
            {
                if (i->excess > 0)
                {
                    *lastActivePtr = i;
                    lastActivePtr = &i->next;
                }
                else i->next = nullptr;
            }
            *lastActivePtr = &nodes[nodeNum];
        }


    template <typename FlowType, typename CostType> 
        inline FlowType SSP<FlowType, CostType>::Augment(Node* start, Node* end)
        {
            FlowType delta = (start->excess < -end->excess) ? start->excess : -end->excess;
            Arc* a;

            for (a=end->parent; a; a=a->sister->head->parent)
            {
                if (delta > a->r_cap) delta = a->r_cap;
            }
            assert(delta > 0);

            end->excess += delta;
            for (a=end->parent; a; a=a->head->parent)
            {
                DecreaseRCap(a, delta);
                a = a->sister;
                IncreaseRCap(a, delta);
            }
            start->excess -= delta;

            return delta;
        }

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::Dijkstra(Node* start)
        {
            assert(start->excess > 0);

            Node* i;
            Node* j;
            Arc* a;
            CostType d;
            Node* permanentNodes;

            std::size_t FLAG0 = ++ counter; // permanently labeled nodes
            std::size_t FLAG1 = ++ counter; // temporarily labeled nodes

            start->parent = nullptr;
            start->flag = FLAG1;
            queue.Reset();
            queue.Add(start, 0);

            permanentNodes = nullptr;

            while ( (i=queue.RemoveMin(d)) )
            {
                assert(i != nullptr);
                if (i->excess < 0)
                {
                    FlowType delta = Augment(start, i);
                    mcf_cost += delta*(d - i->pi + start->pi);
                    for (i=permanentNodes; i; i=i->next_permanent) i->pi += d;
                    break;
                }

                i->pi -= d;
                i->flag = FLAG0;
                i->next_permanent = permanentNodes;
                permanentNodes = i;

                for (a=i->firstNonsaturated; a; a=a->next)
                {
                    j = a->head;
                    if (j->flag == FLAG0) continue;
                    d = a->GetRCost();
                    if (j->flag == FLAG1)
                    {
                        if (d >= queue.GetKey(j)) continue;
                        queue.DecreaseKey(j, d);
                    }
                    else
                    {
                        queue.Add(j, d);
                        j->flag = FLAG1;
                    }
                    j->parent = a;
                }

            }
        }

    template <typename FlowType, typename CostType> 
        inline bool SSP<FlowType, CostType>::node_valid(NodeId i) const
        {
            if(i < 0 || i >= nodeNum) { return false; }
            if(nodes[i].firstSaturated != nullptr) {
                if(nodes[i].firstSaturated < arcs) { return false; }
                if(nodes[i].firstSaturated - arcs >= 2*edgeNum) { return false; }
            }
            if(nodes[i].firstNonsaturated != nullptr) {
                if(nodes[i].firstNonsaturated < arcs) { return false; }
                if(nodes[i].firstNonsaturated - arcs >= 2*edgeNum) { return false; }
            }
            return true;
        }

    template <typename FlowType, typename CostType> 
        inline bool SSP<FlowType, CostType>::arc_valid(Arc* a) const
        {
            if(a < arcs || a >= arcs+2*edgeNum) { return false; }
            if(!node_valid(tail(a-arcs))) { return false; }
            if(!node_valid(head(a-arcs))) { return false; }
            if(a->prev == a) { return false; }
            if(a->next == a) { return false; }
            if(a->sister->sister != a) { return false; }
            if(a->next != nullptr && a->next->prev != a) { return false; }
            if(a->prev != nullptr && a->prev->next != a) { return false; }

            Node* a_tail = a->sister->head;
            if(a_tail->firstSaturated == a && a->prev != nullptr) { return false; }
            if(a_tail->firstNonsaturated == a && a->prev != nullptr) { return false; }

            return true;
        }

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::exchange(Arc* const a, Arc* const b) 
        {
            assert(a >= arcs && a-arcs < 2*edgeNumMax);
            assert(b >= arcs && b-arcs < 2*edgeNumMax);
            assert(arc_valid(a) && arc_valid(b));
            if ( a != b) {
                Arc *sa = a->sister;
                Arc *sb = b->sister;

                Node* a_tail = a->sister->head;
                Node* b_tail = b->sister->head;

                Arc* na = a->next;
                Arc* nb = b->next;

                Arc* pa = a->prev;
                Arc* pb = b->prev;

                if(na) { assert(arc_valid(na)); }
                if(nb) { assert(arc_valid(nb)); }
                if(pa) { assert(arc_valid(pa)); }
                if(pb) { assert(arc_valid(pb)); }

                std::swap(a->head, b->head);
                std::swap(a->r_cap, b->r_cap);
                std::swap(a->cost, b->cost);
                std::swap(capacity[a-arcs], capacity[b-arcs]);

                if ( a != sb) {
                    assert(b != sa);
                    std::swap(a->sister, b->sister);
                    sa->sister = b;
                    sb->sister = a;
                }

                if(a->next == b) {
                    assert(b->prev == a); 
                    auto* b_next = b->next;
                    auto* a_prev = a->prev;
                    b->next = a;
                    a->prev = b;
                    a->next = b_next;
                    b->prev = a_prev;
                    if(b_next) { b_next->prev = a; }
                    if(a_prev) { a_prev->next = b; }
                } else if(b->next == a) {
                    assert(a->prev == b);
                    auto* a_next = a->next;
                    auto* b_prev = b->prev;
                    b->prev = a;
                    a->next = b;
                    a->prev = b_prev;
                    b->next = a_next;
                    if(a_next) { a_next->prev = b; }
                    if(b_prev) { b_prev->next = a; }
                } else {
                    std::swap(a->next, b->next);
                    std::swap(a->prev, b->prev); 
                    if(na != nullptr) { na->prev = b; }
                    if(nb != nullptr) { nb->prev = a; }
                    if(pa != nullptr) { pa->next = b; }
                    if(pb != nullptr) { pb->next = a; } 
                }

                if(a_tail != b_tail) {
                    if(a_tail->firstSaturated == a) { a_tail->firstSaturated = b; }
                    if(a_tail->firstNonsaturated == a) { a_tail->firstNonsaturated = b; }

                    if(b_tail->firstSaturated == b) { b_tail->firstSaturated = a; }
                    if(b_tail->firstNonsaturated == b) { b_tail->firstNonsaturated = a; }
                } else {
                    if(a_tail->firstSaturated == a) { a_tail->firstSaturated = b; }
                    else if(b_tail->firstSaturated == b) { a_tail->firstSaturated = a; }

                    if(a_tail->firstNonsaturated == a) { a_tail->firstNonsaturated = b; }
                    else if(b_tail->firstNonsaturated == b) { a_tail->firstNonsaturated = a; }
                }

                if(na) { assert(arc_valid(na)); }
                if(nb) { assert(arc_valid(nb)); }
                if(pa) { assert(arc_valid(pa)); }
                if(pb) { assert(arc_valid(pb)); }
                assert(arc_valid(a) && arc_valid(b) && arc_valid(sa) && arc_valid(sb));
            }
        }

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::order_inter_nodes()
        {
            for(long e=0; e<2*edgeNum; ++e) { assert(arc_valid(&arcs[e])); }

            std::vector<long> arc_first(nodeNum+1,0);
            std::vector<long> outgoing_arc_index(nodeNum,0);
            for(long e=0; e<2*edgeNum; ++e) {
                arc_first[tail(e)+1]++;
                outgoing_arc_index[tail(e)]++;
            }

            std::partial_sum(arc_first.begin(), arc_first.end(), arc_first.begin());
            std::partial_sum(outgoing_arc_index.begin(), outgoing_arc_index.end(), outgoing_arc_index.begin());

            for(long i=0; i<nodeNum-1; ++i) {

                long last = outgoing_arc_index[i];
                assert(last < 2*edgeNum);

                for ( long arc_num = arc_first[i]; arc_num < last; arc_num ++ ) {
                    long tail_node_id = tail(arc_num);

                    while ( tail_node_id != i ) {
                        long arc_new_num = arc_first[tail_node_id];
                        exchange(&arcs[arc_num], &arcs[arc_new_num]);

                        arc_first[tail_node_id]++;

                        tail_node_id = tail(arc_num);
                    }
                }
            }
        }

    template <typename FlowType, typename CostType> 
        inline void SSP<FlowType, CostType>::order_intra_nodes()
        {
            // sort outgoing arcs of every node by head node id. This assumes that order_inter_nodes has already been called
            std::vector<long> perm(nodeNum);
            std::vector<long> outgoing_arc_begin(nodeNum+1, 0);
            for(EdgeId e=0; e<2*edgeNum; ++e) {
                outgoing_arc_begin[tail(e)+1]++;
            }
            std::partial_sum(outgoing_arc_begin.begin(), outgoing_arc_begin.end(), outgoing_arc_begin.begin());

            for(NodeId i=0; i<nodeNum; ++i) {
                const long no_nodes = outgoing_arc_begin[i+1] - outgoing_arc_begin[i];
                std::iota(perm.begin(), perm.begin() + no_nodes, 0);
                Arc* arc = arcs + outgoing_arc_begin[i];
                Node* n = nodes;

                auto sort_func = [arc,n](auto i, auto j) { 
                    assert(arc[i].sister->head == arc[j].sister->head);
                    return arc[i].head - n < arc[j].head - n; 
                };
                std::sort(perm.begin(), perm.begin()+no_nodes, sort_func);

                // follow cycles in permutation.
                for(long c=0; c<no_nodes; ++c) {
                    long next_idx = perm[c];
                    if(next_idx == c || next_idx < 0) {
                        continue;
                    }
                    long cur_idx = c;
                    while(perm[next_idx] >= 0) {
                        exchange(&arc[cur_idx], &arc[next_idx]); 
                        perm[cur_idx] -= no_nodes; // mark as visited
                        cur_idx = next_idx;
                        next_idx = perm[next_idx];
                    }
                }
                assert(std::is_sorted(arc, arc+no_nodes, [](auto& a, auto& b) { return a.head < b.head; }));
            } 
        }

    template <typename FlowType, typename CostType> 
        inline CostType SSP<FlowType, CostType>::solve()
        {
            assert( 0 == std::accumulate(nodes, nodes+no_nodes(), 0, [](const long int s, const Node& i) { return s + i.excess; }) );
            Node* i;
            Init();
            while ( 1 )
            {
                i = firstActive;
                if (i == &nodes[nodeNum]) break;
                firstActive = i->next;
                i->next = NULL;
                if (i->excess > 0)
                {
                    Dijkstra(i);
                    if (i->excess > 0 && !i->next) 
                    { 
                        i->next = firstActive; 
                        firstActive = i; 
                    }
                }
            }

            assert(TestCosts());
            assert(TestOptimality());

            for(EdgeId e=0; e<2*edgeNum; ++e) { assert(arc_valid(&arcs[e])); }
            return mcf_cost;
        }

    template <typename FlowType, typename CostType> 
        inline CostType SSP<FlowType, CostType>::objective() const
        {
            CostType c = 0.0;
            for(EdgeId a=0; a<2*edgeNum; ++a) {
                c += flow(a)*(arcs[a].cost);
            }
            return c/2.0;
        }

    template <typename FlowType, typename CostType> 
        bool SSP<FlowType, CostType>::TestOptimality() const
        {
            Node* i;
            Arc* a;

            for (i=nodes; i<nodes+nodeNum; i++)
            {
                if (i->excess != 0)
                {
                    return false;
                }
                for (a=i->firstSaturated; a; a=a->next)
                {
                    if (a->r_cap != 0)
                    {
                        return false;
                    }
                }
                for (a=i->firstNonsaturated; a; a=a->next)
                {
                    CostType c = a->GetRCost();
                    if (a->r_cap <= 0 || a->GetRCost() < -1e-5)
                    {
                        return false;
                    }
                }
            }
            return true;
        }

    template <typename FlowType, typename CostType> 
        bool SSP<FlowType, CostType>::TestCosts() const
        {
            CostType _cost = 0;

            for (Arc* a=arcs; a<arcs+2*edgeNum; ++a)
            {
                if(a->r_cap + a->sister->r_cap != capacity[a-arcs] + capacity[a->sister-arcs]) {
                    return false;
                }
            }

            if(std::abs(objective() - mcf_cost)/no_edges() > 1e-8) {
                return false;
            }

            return true;
        }

    template<typename FlowType, typename CostType>
        void SSP<FlowType, CostType>::print_flow() const
        {
            std::cout << "flow:\n";
            for(EdgeId e=0; e<2*edgeNum; ++e) {
                std::cout << tail(e) << " -> " << head(e) << ": flow = " << flow(e) << ", capacity = " << capacity[e] << ", cost = " << arcs[e].cost << "\n";
            }
        }


    // read file in DIMACS format
    template<typename FLOW_TYPE, typename COST_TYPE>
        SSP<FLOW_TYPE,COST_TYPE>* read_dimacs_file(const std::string& filename)
        {
            SSP<FLOW_TYPE,COST_TYPE>* f = nullptr;

            std::ifstream instance;
            instance.open(filename);
            if(!instance.is_open()) {
                throw std::runtime_error("could not open file " + filename);
            }
            std::string line;
            while(std::getline(instance, line))
            {
                // DIMACS format:
                // c arc has <tail> <head> <capacity l.b.> <capacity u.b> <cost>
                if(line.empty()) continue;
                std::istringstream iss(line);
                char id;
                if (!(iss >> id )) { throw std::runtime_error("in file " + filename + ": cannot read line " + line); } 
                switch(id) {
                    case 'c':
                        break;
                    case 'p': 
                        {
                            if(f != nullptr) { throw std::runtime_error("in file " + filename + ": not more than one line beginning with 'p' allowed"); }
                            std::string min;
                            std::size_t n; // number of nodes
                            std::size_t m; // number of arcs
                            if( !(iss >> min >> n >> m)) { throw std::runtime_error("in file " + filename + ": cannot read number of nodes and arcs from line:\n " + line); } 
                            if("min" != min) { throw std::runtime_error("in file " + filename + ": min must come after 'p' in line:\n " + line); } 
                            f = new SSP<FLOW_TYPE,COST_TYPE>(n,m);
                            break;
                        }
                    case 'n': 
                        {
                            std::size_t id, flow;
                            if( !(iss >> id >> flow)) { throw std::runtime_error("in file " + filename + ": cannot read number node id and external flow from line:\n " + line); } 
                            assert(id >= 1);
                            f->add_node_excess(id-1,flow);
                            break;
                        }
                    case 'a': 
                        {
                            long i,j,lower,upper,cost;
                            if( !(iss >> i >> j >> lower >> upper >> cost)) { throw std::runtime_error("in file " + filename + ": cannot read arc information from line:\n " + line); } 
                            assert(i >= 1);
                            assert(j >= 1);
                            f->add_edge(i-1,j-1,lower,upper,cost);
                            break;
                        }
                    default:
                        throw std::runtime_error("in file " + filename + ": unknown line identifier " + std::to_string(id));
                        break;
                }
            }
            if(f == nullptr) { throw std::runtime_error("Exactly one line beginning with 'p' allowed in file " + filename); }
            return f; 
        }

} // namespace MCF

#endif // MCF_SSP_HXX
