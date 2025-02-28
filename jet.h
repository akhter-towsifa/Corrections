#pragma once

#include "correction.h"
#include "corrections.h"
#include "JMESystematicsCalculators.h"

namespace correction {
class JetCorrProvider : public CorrectionsBase<JetCorrProvider> {
public:
     enum class UncSource : int {
        Central = -1,
        JER = 0,
        Total = 1,
        RelativeBal = 2,
        HF = 3,
        BBEC1 = 4,
        EC2 = 5,
        Absolute = 6,
        FlavorQCD = 7,
        BBEC1_year = 8,
        Absolute_year = 9,
        EC2_year = 10,
        HF_year = 11,
        RelativeSample_year = 12
    };

    static const std::string getFullNameUnc(const std::string source_name, const std::string year, bool need_year){
        return need_year ? source_name+year : source_name;
    }
    static const std::string& getScaleStr(UncScale scale)
    {
        static const std::map<UncScale, std::string> names = {
            { UncScale::Down, "down" },
            { UncScale::Central, "nom" },
            { UncScale::Up, "up" },
        };
        return names.at(scale);
    }
    static const int GetScaleIdx(UncScale scale){
        static const std::map<UncScale, int> scale_indexes = {
            {UncScale::Up, 1},
            {UncScale::Down,2},
            {UncScale::Central,0},
        };
        return scale_indexes.at(scale);
    }
    static const int GetJesIdx(UncSource source, UncScale scale){
        int index = static_cast<int>(source) * 2 + GetScaleIdx(scale);
        return std::max(index, 0);
    }
    static const std::map<UncSource,std::tuple<std::string,bool,bool>> getUncMap (){
        static const std::map<UncSource,std::tuple<std::string,bool,bool>> UncMap = {
            {UncSource::Central, {"Central", false,false}},
            {UncSource::JER, {"JER", false,false}},
            {UncSource::FlavorQCD,{"FlavorQCD",true,false}},
            {UncSource::RelativeBal,{"RelativeBal",true,false}},
            {UncSource::HF,{"HF",true,false}},
            {UncSource::BBEC1,{"BBEC1",true,false}},
            {UncSource::EC2,{"EC2",true,false}},
            {UncSource::Absolute,{"Absolute",true,false}},
            {UncSource::Total,{"Total",true,false}},
            {UncSource::BBEC1_year,{"BBEC1_",true,true}},
            {UncSource::Absolute_year,{"Absolute_",true,true}},
            {UncSource::EC2_year,{"EC2_",true,true}},
            {UncSource::HF_year,{"HF_",true,true}},
            {UncSource::RelativeSample_year,{"RelativeSample_",true,true}},
        };
        return UncMap;
    }

    JetCorrProvider(const std::string& ptResolution,const std::string& ptResolutionSF, const std::string& JesTxtFile, const std::string& year)
    {
        jvc_total.setSmearing(ptResolution, ptResolutionSF, false, true, 0.2, 3);
        jvc_total.setAddHEM2018Issue(year=="2018");
        for (auto& [unc_source ,unc_features] : getUncMap()){
            if(! std::get<1>(unc_features) ) continue;
            std::string jes_name = getFullNameUnc(std::get<0>(unc_features), year, std::get<2>(unc_features));
            jvc_total.addJESUncertainty(jes_name,JetCorrectorParameters{JesTxtFile,jes_name});
        }
    }


    std::map<std::pair<UncSource,UncScale>, RVecLV> getShiftedP4(const RVecF& Jet_pt, const RVecF& Jet_eta, const RVecF& Jet_phi,
                    const RVecF& Jet_mass, const RVecF& Jet_rawFactor, const RVecF& Jet_area,
                    const RVecI& Jet_jetId, const float rho, const RVecI& Jet_partonFlavour,
                    std::uint32_t seed, const RVecF& GenJet_pt, const RVecF& GenJet_eta,
                    const RVecF& GenJet_phi, const RVecF& GenJet_mass, int event) const {
        std::map<std::pair<UncSource,UncScale>, RVecLV> all_shifted_p4;
        auto result = jvc_total.produce(Jet_pt, Jet_eta, Jet_phi, Jet_mass, Jet_rawFactor,
                                    Jet_area, Jet_jetId, rho, Jet_partonFlavour, seed,
                                    GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, event);
        std::vector<UncScale> uncScales={UncScale::Central, UncScale::Up, UncScale::Down};
        for (auto& uncScale : uncScales){
            for (auto & [unc_source ,unc_features] : getUncMap()){
                RVecLV shifted_p4(Jet_pt.size());
                if(unc_source != UncSource::Central && uncScale == UncScale::Central) continue;
                if(unc_source == UncSource::Central && uncScale != UncScale::Central) continue;
                int scale_idx = GetJesIdx(unc_source, uncScale);
                for (int jet_idx= 0 ; jet_idx < Jet_pt.size(); ++jet_idx){
                    shifted_p4[jet_idx] = LorentzVectorM(result.pt(scale_idx)[jet_idx], Jet_eta[jet_idx],
                    Jet_phi[jet_idx], result.mass(scale_idx)[jet_idx]);
                }
            all_shifted_p4.insert({{unc_source, uncScale}, shifted_p4});
            }
        }
        return all_shifted_p4;
    }
    RVecF getResolution(const RVecF& Jet_pt, const RVecF& Jet_eta, const float rho) const {
        return jvc_total.getResolution(Jet_pt, Jet_eta, rho);
    }

    //RVecI getVetoMap(const RVecF& Jet_eta, const RVecF& Jet_phi) const{

    //}

private:
    JetVariationsCalculator jvc_total ;
};

    // run3 code starts here
    // main difference: all corrections are retrieved from json file using correctionlib
    class JetCorrectionProvider : public CorrectionsBase<JetCorrectionProvider>
    {
        public:
        enum class UncSource : int
        {
            Central = -1,
            JER = 0,
            Total = 1,
            RelativeBal = 2,
            HF = 3,
            BBEC1 = 4,
            EC2 = 5,
            Absolute = 6,
            FlavorQCD = 7,
            BBEC1_year = 8,
            Absolute_year = 9,
            EC2_year = 10,
            HF_year = 11,
            RelativeSample_year = 12
        };

        // json_file_name - path to json file with corrections
        // e.g. /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2022_Summer2022/jet_jerc.json.gz

        // jecTag_corrName_algoType

        // jec_tag - a string describing when sample was produced and type of sample
        // e.g. Summer22_22Sep2023_V2_MC

        // algo - type of jet algorithm
        // e.g. AK4PFPuppi
        JetCorrectionProvider(std::string const& json_file_name,
                              std::string const& jetsmear_file_name,
                              std::string const& jec_tag,
                              std::string const& jer_tag,
                              std::string const& algo,
                              std::string const& year,
                              bool is_data,
                              bool use_regrouped,
                              bool apply_compound)
        :   corrset_(CorrectionSet::from_file(json_file_name))
        ,   jersmear_corr_(CorrectionSet::from_file(jetsmear_file_name)->at("JERSmear"))
        ,   corr_jer_sf_(corrset_->at(jer_tag + "_ScaleFactor_" + algo))
        ,   corr_jer_res_(corrset_->at(jer_tag + "_PtResolution_" + algo))
        ,   cmpd_corr_(corrset_->compound().at(jec_tag + "_L1L2L3Res_" + algo))
        ,   is_data_(is_data)
        ,   apply_cmpd_(apply_compound)
        {
            // map with uncertainty sources should only be filled for MC
            if (!is_data_)
            {
                auto const& unc_map_ref = use_regrouped ? unc_map_regrouped : unc_map_total;
                for (auto const& [unc_source, unc_name]: unc_map_ref)
                {
                    std::string full_name = jec_tag;
                    full_name += '_';
                    full_name += unc_name;
                    full_name += '_';
                    if (year_dep_map.at(unc_source))
                    {
                        full_name += year;
                        full_name += '_';
                    }
                    full_name += algo;
                    unc_map_[unc_source] = full_name;
                }
            }
        }

        std::map<std::pair<UncSource, UncScale>, RVecLV> getShiftedP4(RVecF Jet_pt, const RVecF& Jet_eta, const RVecF& Jet_phi, RVecF Jet_mass,
                                                                      const RVecF& Jet_rawFactor, const RVecF& Jet_area, const float rho, int event, bool apply_jer,
                                                                      const RVecF& GenJet_pt = {}, const RVecI& Jet_genJetIdx = {}) const
        {
            std::map<std::pair<UncSource, UncScale>, RVecLV> all_shifted_p4;
            std::vector<UncScale> uncScales = { UncScale::Up, UncScale::Down };

            size_t sz = Jet_pt.size();
            std::vector<float> jer_pt_resolutions(sz);
            RVecLV central_p4(sz);
            for (size_t i = 0; i < sz; ++i)
            {
                // uscaling
                if (apply_cmpd_)
                {
                    Jet_pt[i] *= 1.0 - Jet_rawFactor[i];
                    Jet_mass[i] *= 1.0 - Jet_rawFactor[i];
                }

                if (!is_data_ && apply_jer)
                {
                    // extract jer scale factor and resolution
                    float jer_sf = corr_jer_sf_->evaluate({Jet_eta[i], Jet_pt[i], "nom"});
                    float jer_pt_res = corr_jer_res_->evaluate({Jet_eta[i], Jet_pt[i], rho});
                    jer_pt_resolutions[i] = jer_pt_res;

                    int genjet_idx = Jet_genJetIdx[i];
                    float genjet_pt = genjet_idx != -1 ? GenJet_pt[genjet_idx] : -1.0;
                    float jersmear_factor = jersmear_corr_->evaluate({Jet_pt[i], Jet_eta[i], genjet_pt, rho, event, jer_pt_res, jer_sf});

                    // apply jer smearing (only for MC)
                    Jet_pt[i] *= jersmear_factor;
                    Jet_mass[i] *= jersmear_factor;
                }

                // evaluate and apply compound correction
                if (apply_cmpd_)
                {
                    float cmpd_sf = cmpd_corr_->evaluate({Jet_area[i], Jet_eta[i], Jet_pt[i], rho});
                    Jet_pt[i] *= cmpd_sf;
                    Jet_mass[i] *= cmpd_sf;
                }

                central_p4[i] = LorentzVectorM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
            }

            all_shifted_p4.insert({{UncSource::Central, UncScale::Central}, central_p4});

            // apply uncertainties from uncertainty map
            // this part should not be executed for data
            if (!is_data_)
            {
                for (auto const& uncScale: uncScales)
                {
                    for (auto const& [unc_source, unc_name]: unc_map_)
                    {
                        RVecLV shifted_p4(sz);
                        if (unc_source == UncSource::JER)
                        {
                            if (apply_jer)
                            {
                                for (size_t jet_idx = 0; jet_idx < sz; ++jet_idx)
                                {
                                    float sf = 1.0;
                                    sf += static_cast<int>(uncScale)*jer_pt_resolutions[jet_idx];
                                    shifted_p4[jet_idx] = LorentzVectorM(sf*Jet_pt[jet_idx], Jet_eta[jet_idx], Jet_phi[jet_idx], sf*Jet_mass[jet_idx]);
                                }
                                all_shifted_p4.insert({{unc_source, uncScale}, shifted_p4});
                            }
                        }
                        else
                        {
                            for (size_t jet_idx = 0 ; jet_idx < sz; ++jet_idx)
                            {
                                float sf = 1.0;
                                Correction::Ref corr = corrset_->at(unc_name);
                                float unc = corr->evaluate({Jet_eta[jet_idx], Jet_pt[jet_idx]});
                                sf += static_cast<int>(uncScale)*unc;
                                shifted_p4[jet_idx] = LorentzVectorM(sf*Jet_pt[jet_idx], Jet_eta[jet_idx], Jet_phi[jet_idx], sf*Jet_mass[jet_idx]);
                            }
                            all_shifted_p4.insert({{unc_source, uncScale}, shifted_p4});
                        }
                    }
                }
            }
            return all_shifted_p4;
        }

        RVecF GetResolutions(RVecF pt, RVecF mass, RVecF const& raw_factor, RVecF const& eta, float rho) const
        {
            size_t sz = pt.size();
            RVecF res(sz);
            for (size_t i = 0; i < sz; ++i)
            {
                pt[i] *= 1.0 - raw_factor[i];
                mass[i] *= 1.0 - raw_factor[i];
                float jer_sf = corr_jer_sf_->evaluate({eta[i], pt[i], "nom"});
                float jer_pt_res = corr_jer_res_->evaluate({eta[i], pt[i], rho});
                res[i] = jer_pt_res;
            }
            return res;
        }

        private:
        std::map<UncSource, std::string> unc_map_;
        std::unique_ptr<CorrectionSet> corrset_;
        Correction::Ref jersmear_corr_; // aka shared_ptr<Correction const>, sizeof = 8
        Correction::Ref corr_jer_sf_;
        Correction::Ref corr_jer_res_;
        CompoundCorrection::Ref cmpd_corr_;
        bool is_data_;
        bool apply_cmpd_;

        inline static const std::map<UncSource, std::string> unc_map_total = { { UncSource::Total, "Total" },
                                                                               { UncSource::JER, "JER" } };

        inline static const std::map<UncSource, bool> year_dep_map = { { UncSource::Central, false },
                                                                       { UncSource::JER, false },
                                                                       { UncSource::Total, false },
                                                                       { UncSource::RelativeBal, false },
                                                                       { UncSource::HF, false },
                                                                       { UncSource::BBEC1, false },
                                                                       { UncSource::EC2, false },
                                                                       { UncSource::Absolute, false },
                                                                       { UncSource::FlavorQCD, false },
                                                                       { UncSource::BBEC1_year, true },
                                                                       { UncSource::Absolute_year, true },
                                                                       { UncSource::EC2_year, true },
                                                                       { UncSource::HF_year, true },
                                                                       { UncSource::RelativeSample_year, true } };

        inline static const std::map<UncSource, std::string> unc_map_regrouped = { { UncSource::JER, "JER" },
                                                                                   { UncSource::RelativeBal, "Regrouped_RelativeBal" },
                                                                                   { UncSource::HF, "Regrouped_HF" },
                                                                                   { UncSource::BBEC1, "Regrouped_BBEC1" },
                                                                                   { UncSource::EC2, "Regrouped_EC2" },
                                                                                   { UncSource::Absolute, "Regrouped_Absolute" },
                                                                                   { UncSource::FlavorQCD, "Regrouped_FlavorQCD" },
                                                                                   { UncSource::BBEC1_year, "Regrouped_BBEC1" },
                                                                                   { UncSource::Absolute_year, "Regrouped_Absolute" },
                                                                                   { UncSource::EC2_year, "Regrouped_EC2" },
                                                                                   { UncSource::HF_year, "Regrouped_RelativeStatHF" },
                                                                                   { UncSource::RelativeSample_year, "Regrouped_RelativeSample" } };
    };

} // namespace correction