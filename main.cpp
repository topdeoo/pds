#include <fmt/chrono.h>
#include <fmt/format.h>

#include <boost/program_options.hpp>
#include <chrono>
#include <filesystem>
#include <functional>
#include <iostream>
#include <string>

#include "draw_grid.hpp"
#include "graphio.hpp"
#include "gurobi_solve.hpp"
#include "pds.hpp"
#include "pdssolve.hpp"

void printResult(const pds::PdsState& state) {
    size_t active_count = 0, inactive_count = 0, zero_injection_count = 0,
           observed_count = 0;
    auto filter_active = [&state](auto v) -> bool { return state.isActive(v); };
    if (false) {
        auto active = state.graph().vertices() |
                      ranges::views::filter(filter_active) |
                      ranges::to<std::vector>();
        fmt::print("active nodes: {}\n", active);
    }
    for (const auto& v : state.graph().vertices()) {
        switch (state.activeState(v)) {
            case pds::PmuState::Active:
                ++active_count;
                break;
            case pds::PmuState::Inactive:
                ++inactive_count;
                break;
            case pds::PmuState::Blank:
                break;
        }
        zero_injection_count += state.graph()[v].zero_injection;
        observed_count += state.isObserved(v);
    }
    fmt::print(
        "graph (n={}, m={}, #active={}, #inactive={}, #observed={}, "
        "#zero_injection={})\n",
        state.graph().numVertices(), state.graph().numEdges(), active_count,
        inactive_count, observed_count, zero_injection_count);
    bool feasible = state.allObserved();
    fmt::print("solved: {}\n", feasible);
}

void printGraph(const pds::PowerGrid& graph) {
    std::cout << "graph {\n";
    for (auto v : graph.vertices()) {
        std::cout << graph[v].name << "; ";
    }
    std::cout << "\n";
    for (auto e : graph.edges()) {
        std::cout << graph[graph.source(e)].name << " -- "
                  << graph[graph.target(e)].name << ";\n";
    }
    std::cout << "}\n";
}

auto now() { return std::chrono::high_resolution_clock::now(); }

template <typename T>
auto ms(T time) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(time);
}

struct DrawOptions {
    bool drawInput;
    bool drawSolution;
    bool drawReductions;
    bool drawSubproblems;

    bool drawAny() {
        return drawInput || drawSolution || drawReductions || drawSubproblems;
    }
};

auto getModel(const std::string& name) {
    if (name == "gurobi" || name == "jovanovic2") {
        return pds::modelJovanovicExpanded;
    } else if (name == "jovanovic") {
        return pds::modelJovanovic;
    } else if (name == "brimkov") {
        return pds::modelBrimkov;
    } else if (name == "brimkov2") {
        return pds::modelBrimkovExpanded;
    } else if (name == "azami" || name == "azami-brimkov") {
        return pds::modelAzamiBrimkov;
    } else if (name == "domination") {
        return pds::modelDomination;
    } else {
        throw std::invalid_argument("unknown model " + name);
    }
}

int run(int argc, const char** argv) {
    using namespace pds;
    using std::string, std::vector;
    using namespace std::string_literals;
    namespace po = boost::program_options;
    namespace fs = std::filesystem;
    po::options_description desc("options");
    desc.add_options()("help,h", "show this help")(
        "graph,f", po::value<string>(), "input graph")(
        "outdir,o", po::value<string>()->default_value("out"),
        "output directory")("solver,s",
                            po::value<string>()->default_value("gurobi"),
                            "solve method. Can be any of "
                            "[none,greedy,greedy-degree,branching,gurobi,"
                            "brimkov,jovanovic,domination]")(
        "subproblem,u",
        "split problem into subproblems and solve them individually")(
        "print-solve", "print intermediate solve state")(
        "print-state,p", "print solve state after each step")(
        "time-limit,t", po::value<double>()->default_value(600.0),
        "time limit for gurobi in seconds")(
        "reductions,r", "apply reductions before exact solving")(
        "all-zi,z",
        po::value<bool>()->default_value(false)->implicit_value(true, "true"),
        "consider all vertices zero-inection")(
        "draw,d",
        po::value<vector<string>>()
            ->default_value({"none"s}, "none")
            ->implicit_value({"all"s}, "all")
            ->composing(),
        "can be one of [none,all,input,solution,reductions,subproblems]");
    po::positional_options_description pos;
    pos.add("graph", 1);
    po::variables_map vm;
    po::store(
        po::command_line_parser(argc, argv).options(desc).positional(pos).run(),
        vm);
    po::notify(vm);

    if (vm.count("help")) {
        desc.print(std::cout);
        return 1;
    }

    auto outdir = vm["outdir"].as<string>();
    switch (fs::status(outdir).type()) {
        case fs::file_type::none:
        case fs::file_type::not_found:
            fs::create_directories(outdir);
            break;
        case fs::file_type::directory:
        case fs::file_type::symlink:
            break;
        default:
            fmt::print(stderr, "could not create output directory '{}'",
                       outdir);
            return 1;
    }

    DrawOptions drawOptions{};
    auto drawOption = vm["draw"].as<vector<string>>();
    for (auto d : drawOption) {
        if (d == "none") {
            drawOptions.drawSubproblems = drawOptions.drawReductions =
                drawOptions.drawSolution = drawOptions.drawInput = false;
        } else if (d == "all") {
            drawOptions.drawSubproblems = drawOptions.drawReductions =
                drawOptions.drawSolution = drawOptions.drawInput = true;
        } else if (d == "input") {
            drawOptions.drawInput = true;
        } else if (d == "solution") {
            drawOptions.drawSolution = true;
        } else if (d == "reductions") {
            drawOptions.drawReductions = true;
        } else if (d == "subproblems") {
            drawOptions.drawSubproblems = true;
        } else {
            fmt::print(stderr, "invalid draw option {}\n", d);
            return 1;
        }
    }

    string solverName = vm["solver"].as<string>();
    std::function<SolveResult(PdsState&, double)> solve;
    if (solverName == "branching") {
        solve = [](auto& state, double) {
            return solveBranching(state, true,
                                  greedy_strategies::largestDegree);
        };
    } else if (solverName == "greedy") {
        solve = [&vm](auto& state, double) {
            return solveGreedy(state, vm.count("reductions"),
                               greedy_strategies::largestDegree);
        };
    } else if (solverName == "fast-greedy"s) {
        solve = [](auto& state, double) { return fastGreedy(state, true); };
    } else if (solverName == "greedy-degree"s) {
        solve = [&vm](auto& state, double) {
            return solveGreedy(state, vm.count("reductions"),
                               greedy_strategies::largestDegree);
        };
    } else if (solverName == "greedy-median"s) {
        solve = [&vm](auto& state, double) {
            return solveGreedy(state, vm.count("reductions"),
                               greedy_strategies::medianDegree);
        };
    } else if (solverName == "none") {
        solve = [](auto&, double) {
            return SolveResult{{}, {}, SolveState::Other};
        };
    } else {
        try {
            solve = [&vm, model = getModel(solverName)](auto& state,
                                                        double timeout) {
                return solvePowerDominatingSet(state, vm.count("print-solve"),
                                               timeout, noop_v, model);
            };
        } catch (std::invalid_argument& ex) {
            fmt::print(stderr, "{}", ex.what());
            return 2;
        }
    }
    if (!vm.count("graph")) {
        fmt::print(stderr, "no input given\n");
        return 1;
    }

    PdsState state(
        readAutoGraph(vm["graph"].as<string>(), vm["all-zi"].as<bool>()));
    auto input = state;

    if (drawOptions.drawInput) {
        writePds(state.graph(), fmt::format("{}/0_input.pds", outdir));
    }
    auto printState = [&](const PdsState& state) {
        if (vm.count("print-state")) printResult(state);
    };
    fmt::print("input:\n");
    printState(state);

    auto tSolveStart = now();
    size_t counter = 0;

    auto drawCallback = [&](const pds::PdsState& state,
                            const std::string& name) mutable {
        if (drawOptions.drawReductions) {
            writePds(state.graph(), fmt::format("{}/1_red_{:04}_{}.pds", outdir,
                                                counter, name));
            ++counter;
        }
    };

    if (vm.count("reductions")) {
        fmt::print("applying reductions\n");
        exhaustiveReductions(state, true, drawCallback);
        auto tReductions = now();
        printState(state);
        fmt::print("reductions took {}\n", ms(tReductions - tSolveStart));
    }

    vector subproblems = state.subproblems();
    ranges::sort(
        subproblems,
        [](const pds::PdsState& left, const pds::PdsState& right) -> bool {
            return left.graph().numVertices() < right.graph().numVertices();
        });
    if (drawOptions.drawSubproblems) {
        for (size_t i = 0; auto& subproblem : subproblems) {
            if (!subproblem.allObserved()) {
                writePds(subproblem.graph(),
                         fmt::format("{}/comp_{:03}_0unsolved.pds", outdir, i));
                ++i;
            }
        }
    }

    SolveResult result{state.numActive(), state.numActive(),
                       SolveState::Optimal};

    if (vm.count("subproblem")) {
        for (size_t i = 0; auto& subproblem : subproblems) {
            auto tSub = now();
            fmt::print("solving subproblem {}\n", i);
            printState(subproblem);
            size_t initialActive = subproblem.numActive();
            auto subresult = solve(subproblem, vm["time-limit"].as<double>());
            result.state = combineSolveState(result.state, subresult.state);
            state.applySubsolution(subproblem);
            result.lower +=
                std::max(subresult.lower, initialActive) - initialActive;
            result.upper +=
                std::max(subresult.upper, initialActive) - initialActive;
            auto tSubEnd = now();
            fmt::print("solved subproblem {} in {} ({} active)\n", i,
                       ms(tSubEnd - tSub), subproblem.numActive());
            if (drawOptions.drawSubproblems && drawOptions.drawSolution) {
                writePds(subproblem.graph(),
                         fmt::format("{}/comp_{:03}_2solved.pds", outdir, i));
            }
            ++i;
        }
    } else {
        result = solve(state, vm["time-limit"].as<double>());
    }

    if (result.state == SolveState::Infeasible) {
        auto tSolveEnd = now();
        fmt::print("model proved infeasible after {}\n",
                   ms(tSolveEnd - tSolveStart));
        return 1;
    } else {
        for (auto v : state.graph().vertices()) {
            if (state.isActive(v)) {
                input.setActive(v);
            }
        }
        auto tSolveEnd = now();
        if (drawOptions.drawSolution) {
            writePds(state.graph(),
                     fmt::format("{}/2_solved_preprocessed.pds", outdir));
            writePds(input.graph(), fmt::format("{}/3_solved.pds", outdir));
        }
        fmt::print("solved in {}\n", ms(tSolveEnd - tSolveStart));
        printState(state);
        printState(input);

        return 0;
    }
}

int main(int argc, const char** argv) {
    // std::string filename = argv[1];
    // processBoost(filename);
    run(argc, argv);
    return 0;
}
