#pragma once

#include <array>
#include <tuple>
#include <algorithm>
#include <concepts>
#include <mdspan.hpp>
#include <chem/atom.hpp>
#include <chem/expression.hpp>
#include <chem/util.hpp>


namespace chem {

template <std::equality_comparable ID>
struct Species
{
    using type = Species<ID>;
    using id_type = ID;

    constexpr operator id_type() const {
        return static_cast<id_type>(id);
    }

    const id_type id;
    const double mass;
};

template <typename L, typename R>
struct Equation
{
    const L lhs;
    const R rhs;
};

template <typename L, typename R, typename Rate>
struct Reaction
{
    const Equation<L,R> eqn;
    const Rate rate;
};

template <typename ID, typename T>
static constexpr auto operator||(ID lhs, T rhs) {
    return Species<ID> {lhs, atom::atomic_mass(rhs)};
}

template <typename L, typename R>
static constexpr auto operator>=(L lhs, R rhs) {
    return Equation {lhs, rhs};
}

template <typename L, typename R, typename Rate>
static constexpr auto operator||(Equation<L,R> lhs, Rate rhs) {
    return Reaction {lhs, rhs};
}

template <typename S, size_t N>
struct VariableSpecies : std::array<S, N> 
{ };

template <typename T0, typename... TN>
VariableSpecies(T0, TN...) -> VariableSpecies<typename T0::type, 1+sizeof...(TN)>;

template <typename S, size_t N>
struct FixedSpecies : std::array<S, N> 
{ };

template <typename T0, typename... TN>
FixedSpecies(T0, TN...) -> FixedSpecies<typename T0::type, 1+sizeof...(TN)>;


template <VariableSpecies Var, FixedSpecies Fix, auto... React>
class Mechanism
{
public:

    using species_type = decltype(Var)::value_type;
    using species_id_type = species_type::id_type;

    static constexpr size_t nvar = Var.size();
    static constexpr size_t nfix = Fix.size();
    static constexpr size_t nspc = nvar + nfix;
    static constexpr size_t nrct = sizeof...(React);

    template <typename... Args>
    static constexpr auto rates(Args... args) {
        return std::array<double, nrct> { calc_rate(React.rate, args...)... };
    }

    static constexpr auto lhs_stoich() {
        std::array<double, nrct*nspc> data;
        auto mat = std::experimental::mdspan(data.data(), nrct, nspc);
        data.fill(0);
        constexpr std::tuple react {React...};
        for_constexpr<0, nrct>([&](auto i) {
            constexpr auto terms = equation_terms(std::get<i>(react).eqn.lhs);
            for_constexpr<0, terms.size()>([&](auto n) {
                const size_t j = spc_num<terms[n].first>();
                mat[i,j] += terms[n].second;
            });
        });
        return data;
    }

    static constexpr auto rhs_stoich() {
        std::array<double, nrct*nspc> data;
        auto mat = std::experimental::mdspan(data.data(), nrct, nspc);
        data.fill(0);
        constexpr std::tuple react {React...};
        for_constexpr<0, nrct>([&](auto i) {
            constexpr auto terms = equation_terms(std::get<i>(react).eqn.rhs);
            for_constexpr<0, terms.size()>([&](auto n) {
                const size_t j = spc_num<terms[n].first>();
                mat[i,j] += terms[n].second;
            });
        });
        return data;
    }

    static constexpr auto agg_stoich() {
        std::array<double, nrct*nspc> data;
        auto mat = std::experimental::mdspan(data.data(), nrct, nspc);
        data.fill(0);
        constexpr std::tuple react {React...};
        for_constexpr<0, nrct>([&](auto i) {
            constexpr auto eqn = std::get<i>(react).eqn;
            constexpr auto lhs_terms = equation_terms(eqn.lhs);
            for_constexpr<0, lhs_terms.size()>([&](auto n) {
                constexpr auto id = lhs_terms[n].first;
                constexpr auto coef = lhs_terms[n].second;
                const size_t j = spc_num<id>();
                if constexpr (is_var_spc<id>()) {
                    mat[i,j] -= coef;
                }
            });
            constexpr auto rhs_terms = equation_terms(eqn.rhs);
            for_constexpr<0, rhs_terms.size()>([&](auto n) {
                constexpr auto id = rhs_terms[n].first;
                constexpr auto coef = rhs_terms[n].second;
                const size_t j = spc_num<id>();
                if constexpr (is_var_spc<id>()) {
                    mat[i,j] += coef;
                }
            });
        });
        return data;
    }

    static constexpr auto lhs_stoich_rank() {
        
    }

private:

    using term_type = std::pair<species_id_type, double>;

    static constexpr auto equation_terms(const species_id_type& expr) {
        return std::array {term_type{expr, 1.0}};
    }

    template <Arithmetic A>
    static constexpr auto equation_terms(const Product<species_id_type, A>& expr) {
        return std::array {term_type{expr.lhs, expr.rhs}};
    }

    template <Arithmetic A>
    static constexpr auto equation_terms(const Product<A, species_id_type>& expr) {
        return std::array {term_type{expr.rhs, expr.lhs}};
    }
    
    template <Arithmetic A, typename RHS>
    static constexpr auto equation_terms(const Product<A, RHS>& expr) {
        double coef = expr.lhs;
        auto rhs_terms = equation_terms(expr.rhs);
        std::array<term_type, rhs_terms.size()> terms;
        for (size_t i=0; i<terms.size(); ++i) {
            terms[i] = term_type{rhs_terms[i].first, coef*rhs_terms[i].second};
        }
        return terms;
    }

    template <typename LHS, Arithmetic A>
    static constexpr auto equation_terms(const Product<LHS, A>& expr) {
        auto lhs_terms = equation_terms(expr.lhs);
        double coef = expr.rhs;
        std::array<term_type, lhs_terms.size()> terms;
        for (size_t i=0; i<terms.size(); ++i) {
            terms[i] = term_type{lhs_terms[i].first, coef*lhs_terms[i].second};
        }
        return terms;
    }

    template <typename LHS, typename RHS>
    static constexpr auto equation_terms(const Sum<LHS, RHS>& sum) {
        auto lhs_terms = equation_terms(sum.lhs);
        auto rhs_terms = equation_terms(sum.rhs);
        std::array<term_type, lhs_terms.size()+rhs_terms.size()> terms;
        std::copy(rhs_terms.begin(), rhs_terms.end(), std::copy(lhs_terms.begin(), lhs_terms.end(), terms.begin()));
        return terms;
    }

    template <Arithmetic Rate, typename... Args>
    static constexpr double calc_rate(const Rate& rconst, Args&&...) {
        return rconst;
    }

    template <typename Rate, typename... Args>
    requires std::is_invocable_v<Rate, Args...>
    static constexpr double calc_rate(const Rate& rfun, Args&&... args) {
        return rfun(args...);
    }

    template <typename Rate, typename... Args>
    static constexpr double calc_rate(const Rate&, Args&&...) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    template <species_id_type ID>
    static constexpr size_t spc_num() {
        for (size_t i=0; i<spc_index.size(); ++i) {
            if (spc_index[i] == ID) {
                return i;
            }
        }
        return spc_index.size();
    }

    template <species_id_type ID>
    static constexpr bool is_var_spc() {
        for (size_t i=0; i<Var.size(); ++i) {
            if (static_cast<species_id_type>(Var[i]) == ID) {
                return true;
            }
        }
        return false;
    }
    
    static constexpr auto declared_order() {
        std::array<species_id_type, nspc> order;
        for (size_t i=0; i<Var.size(); ++i) {
            order[i] = static_cast<species_id_type>(Var[i]);
        }
        for (size_t i=0; i<Fix.size(); ++i) {
            order[nvar+i] = static_cast<species_id_type>(Fix[i]);
        }
        return order;
    }

    static constexpr auto spc_index = declared_order();
};

} // namespace Chem