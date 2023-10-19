//
// Created by max on 19.07.22.
//

#ifndef PDS_PDS_HPP
#define PDS_PDS_HPP

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <fstream>
#include <iostream>
#include <mpgraphs/graph.hpp>
#include <mpgraphs/vecmap.hpp>
#include <mpgraphs/vecset.hpp>
#include <optional>
#include <range/v3/all.hpp>

#include "map.hpp"
#include "utility.hpp"

namespace pds {

enum class SolveState { Optimal, Timeout, Infeasible, Heuristic, Other };

struct SolveResult {
    size_t lower;
    size_t upper;
    SolveState state;
};

SolveState combineSolveState(SolveState first, SolveState second);

enum class PmuState : signed char {
    Blank = 0,
    Active = 1,
    Inactive = -1
};  // TAG Active means in pre-selected and solution, Blank means undecided,
    // Inactive means excluded

struct Bus {
    std::string name;
    long id;
    bool zero_injection;  // TAG unknown use
    PmuState pmu;
};

#ifdef USE_HASHMAP
using PowerGrid = mpgraphs::MapGraph<Bus, mpgraphs::Empty,
                                     mpgraphs::EdgeDirection::Undirected>;
using ObservationGraph =
    mpgraphs::MapGraph<mpgraphs::Empty, mpgraphs::Empty,
                       mpgraphs::EdgeDirection::Bidirectional>;
using VertexSet = mpgraphs::set<PowerGrid::VertexDescriptor>;
template <typename T>
using VertexMap = mpgraphs::map<PowerGrid::VertexDescriptor, T>;
#else
using Timestamp = std::uint8_t;
using PowerGrid = mpgraphs::VecGraph<Bus, mpgraphs::EdgeDirection::Undirected,
                                     true, Timestamp, std::uint32_t>;
using ObservationGraph =
    mpgraphs::VecGraph<mpgraphs::Empty, mpgraphs::EdgeDirection::Bidirectional,
                       true, Timestamp>;
using VertexSet = mpgraphs::VecSet<PowerGrid::VertexDescriptor, Timestamp>;
template <typename T>
using VertexMap = mpgraphs::VecMap<PowerGrid::VertexDescriptor, T, Timestamp>;
#endif
using VertexList = std::vector<PowerGrid::VertexDescriptor>;

template <class T>
size_t intersectionSize(const set<T>& first, const set<T>& second) {
    size_t count = 0;
    for (auto x : first) {
        if (second.contains(x)) {
            ++count;
        }
    }
    return count;
}

template <class T>
bool isSuperset(const set<T>& container, const set<T>& other) {
    if (other.size() > container.size()) return false;
    for (auto x : other) {
        if (!container.contains(x)) return false;
    }
    return true;
}

class PdsState {
   public:
    using Vertex = PowerGrid::VertexDescriptor;

   private:
    VertexMap<ssize_t> m_unobserved_degree;  // 顶点 v 的邻居中未被观察的顶点数
    VertexSet m_seen;      // 已被访问的顶点集合（dfs用）
    size_t m_numActive;    // Active状态的顶点数
    size_t m_numInactive;  // Inactive状态的顶点数
    PowerGrid m_graph;     //  G = <V, E>
    ObservationGraph
        m_dependencies;  // observed关系依赖图（包括propagating所产生的依赖关系）
    // mpgraphs::MapGraph<mpgraphs::Empty, mpgraphs::Empty,
    // mpgraphs::EdgeDirection::Bidirectional> m_dependencies;

    void propagate(std::vector<Vertex>& queue);

    bool observe(Vertex vertex, Vertex origin);
    bool observeOne(Vertex vertex, Vertex origin, std::vector<Vertex>& queue);

    inline bool disableLowDegreeRecursive(PowerGrid::VertexDescriptor start,
                                          VertexSet& seen);

   public:
    PdsState();
    explicit PdsState(PowerGrid&& graph);

    explicit PdsState(const PowerGrid& graph);

    PdsState(const PdsState&) = default;
    PdsState(PdsState&&) = default;
    PdsState& operator=(const PdsState&) = default;
    PdsState& operator=(PdsState&&) = default;

    void addEdge(Vertex source, Vertex target);
    void removeVertex(Vertex v);

    [[nodiscard]] inline bool isZeroInjection(Vertex vertex) const {
        return m_graph[vertex].zero_injection;
    }

    [[nodiscard]] inline PmuState activeState(Vertex vertex) const {
        return m_graph[vertex].pmu;
    }

    [[nodiscard]] inline bool isActive(Vertex vertex) const {
        return activeState(vertex) == PmuState::Active;
    }

    [[nodiscard]] inline bool isBlank(Vertex vertex) const {
        return activeState(vertex) == PmuState::Blank;
    }

    [[nodiscard]] inline bool isInactive(Vertex vertex) const {
        return activeState(vertex) == PmuState::Inactive;
    }

    [[nodiscard]] inline bool isObserved(Vertex vertex) const {
        return m_dependencies.hasVertex(vertex);
    }

    [[nodiscard]] inline bool isObservingEdge(Vertex source,
                                              Vertex target) const {
        return m_dependencies.hasEdge(source, target);
    }

    [[nodiscard]] inline size_t unobservedDegree(Vertex v) const {
        assert(m_unobserved_degree.at(v) ==
               ranges::distance(m_graph.neighbors(v) |
                                ranges::views::filter([this](auto v) {
                                    return !isObserved(v);
                                })));
        return m_unobserved_degree.at(v);
    }

    bool setActive(Vertex vertex);

    bool unsetActive(Vertex vertex);

    bool setInactive(Vertex vertex);

    bool setBlank(Vertex vertex);

    [[nodiscard]] inline bool allObserved() const {
        return numObserved() == graph().numVertices();
    }

    [[nodiscard]] inline size_t numObserved() const {
        return m_dependencies.numVertices();
    }

    [[nodiscard]] inline size_t numUnobserved() const {
        return graph().numVertices() - numObserved();
    }

    [[nodiscard]] inline size_t numActive() const {
        assert(ranges::distance(graph().vertices() |
                                ranges::views::filter([this](auto v) {
                                    return isActive(v);
                                })) == (ssize_t)m_numActive);
        return m_numActive;
    }

    [[nodiscard]] inline size_t numInactive() const {
        assert(ranges::distance(graph().vertices() |
                                ranges::views::filter([this](auto v) {
                                    return isInactive(v);
                                })) == (ssize_t)m_numInactive);
        return m_numInactive;
    }

    [[nodiscard]] inline size_t numBlank() const {
        return graph().numVertices() - numActive() - numInactive();
    }

    [[nodiscard]] inline size_t numZeroInjection() const {
        return ranges::distance(graph().vertices() |
                                ranges::views::filter([this](auto v) {
                                    return isZeroInjection(v);
                                }));
    }

    [[nodiscard]] inline const PowerGrid& graph() const { return m_graph; }
    [[nodiscard]] inline const auto& observationGraph() const {
        return m_dependencies;
    }

    inline PowerGrid&& moveGraph() { return std::move(m_graph); }

    bool collapseLeaves();

    bool disableLowDegree();

    bool reduceObservedNonZi();

    bool collapseDegreeTwo();

    bool disableObservationNeighborhood();

    bool activateNecessaryNodes();

    bool collapseObservedEdges();

    [[nodiscard]] std::vector<PdsState> subproblems() const;

    void applySubsolution(const PdsState& other);
};

}  // namespace pds

#endif  // PDS_PDS_HPP
