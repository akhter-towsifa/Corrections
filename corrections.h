#pragma once

namespace correction {

using LorentzVectorXYZ = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>;
using LorentzVectorM = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;
using LorentzVectorE = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double>>;
using RVecI = ROOT::VecOps::RVec<int>;
using RVecS = ROOT::VecOps::RVec<size_t>;
using RVecShort = ROOT::VecOps::RVec<short>;
using RVecUC = ROOT::VecOps::RVec<UChar_t>;
using RVecF = ROOT::VecOps::RVec<float>;
using RVecB = ROOT::VecOps::RVec<bool>;
using RVecVecI = ROOT::VecOps::RVec<RVecI>;
using RVecLV = ROOT::VecOps::RVec<LorentzVectorM>;

enum class UncScale : int {
    Down = -1,
    Central = 0,
    Up = +1,
};


template<typename ...Args>
void print_args(Args&&... args) noexcept
{
   ((std::cerr << "\t" << std::forward<Args>(args) << "\n"), ...);
}

template <typename CorrectionClass>
class CorrectionsBase {
public:
    template<typename ...Args>
    static void Initialize(Args&&... args)
    {
        try
        {
            _getGlobal() = std::make_unique<CorrectionClass>(args...);
        }
        catch (std::exception& e)
        {
            std::cerr << "Erorr while initializing " << typeid(CorrectionClass).name() << " with arguments:" << "\n";
            print_args(std::forward<Args>(args)...);
            std::cerr << "exception category: " << e.what() << "\n";
            throw;
        }
        catch (...)
        {
            std::cerr << "Erorr while initializing " << typeid(CorrectionClass).name() << " with arguments:" << "\n";
            print_args(std::forward<Args>(args)...);
            std::cerr << "exception category: unknown\n";
            throw;
        }
    }

    static const CorrectionClass& getGlobal()
    {
        const auto& corr = _getGlobal();
        if(!corr)
        {
            std::cerr << typeid(CorrectionClass).name() << " used before initialization";
            throw std::runtime_error("Class not initialized.");
        }
        return *corr;
    }

private:
    static std::unique_ptr<CorrectionClass>& _getGlobal()
    {
        static std::unique_ptr<CorrectionClass> corr;
        return corr;
    }
};
} // end of namespace correction