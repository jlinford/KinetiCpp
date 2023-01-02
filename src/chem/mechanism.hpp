#pragma once

#include <vector>
#include <functional>


namespace chem {


enum struct Atom : size_t
{
    // Atomic number 0 has zero mass, i.e. ignore mass balance
    ZMB = 0,
    // Atomic numbers
    H,
    He,
    Li,
    Be,
    B,
    C,
    N,
    O,
    F,
    Ne
};


struct Species
{
    // Atomic composition could be useful for mass balance checks
    using AtomList = std::vector<Atom>;
    const AtomList atoms;
    
    // Species number is just a column index in the stoichiometric matrix
    size_t number;

    bool operator==(const Species & rhs) {
        return number == rhs.number;
    }
};


struct Term
{
    Term(const Species & _spc) : 
        coef(1.0), spc(_spc) 
    {}

    Term(const double _coef, const Species & _spc) :
        coef(_coef), spc(_spc)
    {}

    const double coef;
    const Species & spc;
};

struct Reaction
{
    using Expression = std::vector<Term>;
    using RateFunction = std::function<double(double)>;

    const Expression lhs;
    const Expression rhs;
    const RateFunction rate;
};



class Mechanism
{
    using SpeciesList = std::vector<std::reference_wrapper<Species>>;
    using ReactionList = std::vector<Reaction>;
    
    const SpeciesList var_spc_;
    const SpeciesList fix_spc_;
    const ReactionList react_;

    void IndexSpecies() 
    {
        size_t num = 0;
        for (Species & spc : var_spc_) {
            spc.number = num;
            ++num;
        }
        for (Species & spc : fix_spc_) {
            spc.number = num;
            ++num;
        }
    }

public:

    Mechanism(
        const SpeciesList && var_spc, 
        const SpeciesList && fix_spc,
        const ReactionList && react) :
            var_spc_(var_spc),
            fix_spc_(fix_spc),
            react_(react)
        {
            IndexSpecies();
        }

    const size_t nvar() const {
        return var_spc_.size();
    }

    const size_t nfix() const {
        return fix_spc_.size();
    }

    const size_t nspec() const {
        return nvar() + nfix();
    }

    const size_t nreact() const {
        return react_.size();
    }
}; // Mechanism

} // namespace Chem