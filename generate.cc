#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <lemon/core.h>
#include <lemon/connectivity.h>
#include <lemon/euler.h>
#include <lemon/maps.h>
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>
#include <lemon/cost_scaling.h>
#include <lemon/capacity_scaling.h>
#include <lemon/cycle_canceling.h>

using namespace std;

typedef lemon::ListDigraph G;

struct Word {
    G::Node suffix, prefix;
    G::Node tour_node;
    G::Arc extra, final, initial;
    bool redundant;
};

struct Edge {
    unordered_map<string, Word>::iterator w;
    G::Arc arc;
};

struct Affix {
    vector<Edge> suffix, prefix;
    G::Node node;
    G::Node tour_node;
};

template<class MCF>
bool solve(const G &net, const G::ArcMap<int> &lowerMap, const G::ArcMap<int> &upperMap, const G::ArcMap<int> &costMap, const G::NodeMap<int> &supplyMap, int &totalCost, G::ArcMap<int> &flowMap)
{
    MCF mcf(net);
    if (mcf.lowerMap(lowerMap).upperMap(upperMap).costMap(costMap).supplyMap(supplyMap).run() != mcf.OPTIMAL)
        return false;
    totalCost = mcf.totalCost();
    mcf.flowMap(flowMap);
    return true;
}

G::Node &makeNode(G &g, G::Node &n)
{
    if (n == lemon::INVALID)
        n = g.addNode();
    return n;
}

int main(int argc, char **argv)
{
    clog << "Reading dictionary from stdin" << endl;
    unordered_map<string, Affix> affixes;
    unordered_map<string, Word> words;
    unordered_set<string> subwords;
    G net, tour;
    G::ArcMap<int> lowerMap(net), upperMap(net), costMap(net);
    G::NodeMap<int> supplyMap(net);
    string new_word;
    while (getline(cin, new_word)) {
        bool redundant = subwords.find(new_word) != subwords.end();
        for (auto i = new_word.begin(); i != new_word.end(); ++i) {
            for (auto j = new_word.end(); j != i; --j) {
                string s(i, j);
                auto w = words.find(s);
                if (w != words.end())
                    w->second.redundant = true;
                subwords.insert(s);
            }
        }
        words.emplace(new_word, Word {.suffix = lemon::INVALID, .prefix = lemon::INVALID, .tour_node = lemon::INVALID, .extra = lemon::INVALID, .final = lemon::INVALID, .initial = lemon::INVALID, .redundant = redundant});
    }
    G::Node initial = net.addNode(), final = net.addNode();
    supplyMap.set(initial, 1);
    supplyMap.set(final, -1);
    for (auto w = words.begin(); w != words.end(); ++w) {
        w->second.suffix = net.addNode();
        supplyMap.set(w->second.suffix, w->second.redundant ? 0 : 1);
        w->second.prefix = net.addNode();
        supplyMap.set(w->second.prefix, w->second.redundant ? 0 : -1);
        w->second.extra = net.addArc(w->second.prefix, w->second.suffix);
        lowerMap.set(w->second.extra, 0);
        upperMap.set(w->second.extra, numeric_limits<int>::max());
        costMap.set(w->second.extra, 0);
        for (auto i = next(w->first.begin()); i != w->first.end(); ++i) {
            affixes.emplace(string(w->first.begin(), i), Affix()).first->second.prefix.push_back(Edge {w});
            affixes.emplace(string(i, w->first.end()), Affix()).first->second.suffix.push_back(Edge {w});
        }
        if (!w->second.redundant) {
            w->second.final = net.addArc(w->second.suffix, final);
            lowerMap.set(w->second.final, 0);
            upperMap.set(w->second.final, 1);
            costMap.set(w->second.final, 0);
            w->second.initial = net.addArc(initial, w->second.prefix);
            lowerMap.set(w->second.initial, 0);
            upperMap.set(w->second.initial, 1);
            costMap.set(w->second.initial, w->first.length());
        }
    }
    for (auto a = affixes.begin(); a != affixes.end();) {
        if (a->second.suffix.empty() || a->second.prefix.empty() ||
            (a->second.suffix.size() == 1 && a->second.prefix.size() == 1 &&
             a->second.suffix.begin()->w == a->second.prefix.begin()->w)) {
            affixes.erase(a++);
        } else {
            a->second.node = net.addNode();
            supplyMap.set(a->second.node, 0);
            for (Edge &e : a->second.suffix) {
                e.arc = net.addArc(e.w->second.suffix, a->second.node);
                lowerMap.set(e.arc, 0);
                upperMap.set(e.arc, numeric_limits<int>::max());
                costMap.set(e.arc, 0);
            }
            for (Edge &e : a->second.prefix) {
                e.arc = net.addArc(a->second.node, e.w->second.prefix);
                lowerMap.set(e.arc, 0);
                upperMap.set(e.arc, numeric_limits<int>::max());
                costMap.set(e.arc, e.w->first.length() - a->first.length());
            }
            a->second.tour_node = lemon::INVALID;
            ++a;
        }
    }

    clog << "Read " << words.size() << " words and found " << affixes.size() << " affixes; ";
    clog << "created network with " << countNodes(net) << " nodes and " << countArcs(net) << " arcs" << endl;

    int totalCost;
    G::ArcMap<int> flowMap(net);
    bool solved;
    if (argc > 1 && string(argv[1]) == "--net") {
        clog << "Using network simplex algorithm" << endl;
        solved = solve<lemon::NetworkSimplex<G>>(net, lowerMap, upperMap, costMap, supplyMap, totalCost, flowMap);
    } else if (argc > 1 && string(argv[1]) == "--cap") {
        clog << "Using capacity scaling algorithm" << endl;
        solved = solve<lemon::CapacityScaling<G>>(net, lowerMap, upperMap, costMap, supplyMap, totalCost, flowMap);
    } else if (argc > 1 && string(argv[1]) == "--cycle") {
        clog << "Using cycle canceling algorithm" << endl;
        solved = solve<lemon::CycleCanceling<G>>(net, lowerMap, upperMap, costMap, supplyMap, totalCost, flowMap);
    } else if ((argc > 1 && string(argv[1]) == "--cost") || true) {
        clog << "Using cost scaling algorithm" << endl;
        solved = solve<lemon::CostScaling<G>>(net, lowerMap, upperMap, costMap, supplyMap, totalCost, flowMap);
    }

    if (!solved) {
        clog << "error: no solution found" << endl;
        return 1;
    }
    clog << "Lower bound: " << totalCost << endl;

    G::ArcMap<string> arcLabel(tour);
    G::Node empty = tour.addNode();
    for (auto &w : words) {
        if (flowMap[w.second.final])
            arcLabel.set(tour.addArc(makeNode(tour, w.second.tour_node), empty), "");
        if (flowMap[w.second.initial])
            arcLabel.set(tour.addArc(empty, makeNode(tour, w.second.tour_node)), w.first);
    }
    for (auto &a : affixes) {
        for (Edge &e : a.second.suffix) {
            for (int i = 0; i < flowMap[e.arc]; i++)
                arcLabel.set(tour.addArc(makeNode(tour, e.w->second.tour_node), makeNode(tour, a.second.tour_node)), "");
        }
        for (Edge &e : a.second.prefix) {
            for (int i = 0; i < flowMap[e.arc]; i++)
                arcLabel.set(tour.addArc(makeNode(tour, a.second.tour_node), makeNode(tour, e.w->second.tour_node)), e.w->first.substr(a.first.length()));
        }
    }

    clog << "Created tour graph with " << countNodes(tour) << " nodes and " << countArcs(tour) << " arcs" << endl;

    if (!lemon::eulerian(tour)) {
        clog << "error: tour graph is not Eulerian" << endl;
        return 1;
    }

    for (lemon::DiEulerIt<G> e(tour, empty); e != lemon::INVALID; ++e)
        cout << arcLabel[e];

    return 0;
}
