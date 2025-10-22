# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Julia 1.5.1
#     language: julia
#     name: julia-1.5
# ---

# # Expr
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

    abstract type SelectionMethod end
    abstract type CrossoverMethod end
    abstract type MutationMethod end

    struct TruncationSelection <: SelectionMethod
        k::Int
    end

    struct RouletteWheelSelection <: SelectionMethod end

        using ExprRules
        using TreeView
        using Distributions
        using Random

        struct TreeCrossover <: CrossoverMethod
            grammar::Grammar
            max_depth::Int
        end
        struct TreeMutation <: MutationMethod
            grammar::Grammar
            p::Float64
        end

        g = let

            grammar = @grammar begin
                R = |(1:9)
                R = R + R
                R = R - R
                R = R / R
                R = R * R
            end

            function select(T::TruncationSelection, y)
                p = sortperm(y)
                return [p[rand(1:T.k, 2)] for i in y]
            end
            function select(::RouletteWheelSelection, y)
                y = maximum(y) - y
                cat = Categorical(normalize(y, 1))
                return [rand(cat, 2) for i in y]
            end

            function crossover(C::TreeCrossover, a, b)
                child = deepcopy(a)
                crosspoint = sample(b)
                typ = return_type(C.grammar, crosspoint.ind)
                d_subtree = depth(crosspoint)
                d_max = C.max_depth + 1 - d_subtree
                if d_max > 0 && contains_returntype(child, C.grammar, typ, d_max)
                    loc = sample(NodeLoc, child, typ, C.grammar, d_max)
                    insert!(child, loc, deepcopy(crosspoint))
                end
                child
            end

            function mutate(M::TreeMutation, a)
                child = deepcopy(a)
                if rand() < M.p
                    # println("mutate!")
                    loc = sample(NodeLoc, child)
                    # @show loc
                    typ = return_type(M.grammar, get(child, loc).ind)
                    # @show loc
                    subtree = rand(RuleNode, M.grammar, typ)
                    # @show subtree
                    child = insert!(child, loc, subtree)
                    # @show child
                end
                return child
            end

            f = (node) -> begin
                value = Core.eval(node, grammar)
                if isinf(value) || isnan(value)
                    return Inf
                end
                Δ = abs(value - π)
                return log(Δ) + length(node)/1e3
            end

            function genetic_algorithm(f, population, k_max, S, C, M)
                for k in 1 : k_max
                    parents = select(S, f.(population))
                    children = [crossover(C,population[p[1]],population[p[2]]) for p in parents]
                    population = [mutate(M, c) for c in children]
                end
                population[argmin(f.(population))]
            end

            Random.seed!(1)
            m = 1000
            population = [rand(RuleNode, grammar, :R) for i in 1:m]
            k_max = 20
            best_tree = genetic_algorithm(f, population, k_max,
                TruncationSelection(50),
                TreeCrossover(grammar, 10),
                TreeMutation(grammar, 0.25))
            expr = get_executable(best_tree, grammar)
            tree = walk_tree(expr)
            picture = TreeView.tikz_representation(tree)
            picture.options = "scale = 0.6"
            picture
        end

        global cur_plot = g
