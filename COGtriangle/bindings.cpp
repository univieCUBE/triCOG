#include "cogmaker.h"

#include <vector>
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>



namespace py = pybind11;

// Helper to convert Python dict to map
std::map<int, std::vector<BlastHitData>> convert_blast_hits(py::dict blast_dict)
{
    std::map<int, std::vector<BlastHitData>> result;
    
    // Iterate through dict
    for (const auto& item : blast_dict) {
        int query = item.first.cast<int>();
        // Get list of dicts
        py::list hits = item.second.cast<py::list>();
        
        std::vector<BlastHitData> hit_vec;
        // Iterate through dicts
        for (auto hit : hits) {
            py::dict h = hit.cast<py::dict>();
            BlastHitData bhd;
            bhd.query = h["query"].cast<int>();
            bhd.subject = h["subject"].cast<int>();
            bhd.qStart = h["qStart"].cast<double>();
            bhd.qEnd = h["qEnd"].cast<double>();
            bhd.sStart = h["sStart"].cast<double>();
            bhd.sEnd = h["sEnd"].cast<double>();
            bhd.evalue = h["evalue"].cast<double>();
            bhd.score = h["score"].cast<double>();
            // Add struct to hit_vec
            hit_vec.push_back(bhd);
        }
        result[query] = hit_vec;
    }
    return result;
}

// Helper to convert Python dict to std::map<int, std::set<int>>
    std::map<int, std::set<int>> convert_query2subject(py::object obj) {
        if (obj.is_none()) {
            return {};
        }
        
        py::dict query2subject_dict = obj.cast<py::dict>();
        std::map<int, std::set<int>> result;
        
        for (const auto& item : query2subject_dict) {
            int query_id = item.first.cast<int>();
            py::set subject_set = item.second.cast<py::set>();
            
            std::set<int> cpp_set;
            for (auto subj : subject_set) {
                cpp_set.insert(subj.cast<int>());
            }
            result[query_id] = cpp_set;
        }
        return result;
    }


// Create python module
PYBIND11_MODULE(cogmaker, m) {
    // Setup logger
    m.def("set_log_level", [](const int level) {
        auto logger = spdlog::stdout_color_mt("cogmaker");
        
        if (!logger){
            std::cout << "Logger is not set" << std::endl;
            return;
        } 
        
        // Add console logger
        auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        
        if (level == 10) {
            // Add file logger for debug level
            auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("debug.log", false);

            // Set sinks
            logger->sinks() = {console_sink, file_sink};
            logger->set_level(spdlog::level::debug);
            
            // Avoid buffering
            logger->flush_on(spdlog::level::debug); 
        } 
        else {
            // Only console logger for higher levels
            logger->sinks() = {console_sink};

            if (level == 20) logger->set_level(spdlog::level::info);
            else if (level == 30) logger->set_level(spdlog::level::warn);
            else logger->set_level(spdlog::level::err);
        }
    });
    
    // Expose Lse struct
    py::class_<Lse>(m, "Lse")
        .def(py::init<>());

    // Expose BlastHitData struct
    py::class_<BlastHitData>(m, "BlastHitData")
        .def(py::init<>())
        .def_readwrite("query", &BlastHitData::query)
        .def_readwrite("subject", &BlastHitData::subject)
        .def_readwrite("qStart", &BlastHitData::qStart)
        .def_readwrite("qEnd", &BlastHitData::qEnd)
        .def_readwrite("sStart", &BlastHitData::sStart)
        .def_readwrite("sEnd", &BlastHitData::sEnd)
        .def_readwrite("evalue", &BlastHitData::evalue)
        .def_readwrite("score", &BlastHitData::score);
    
    // Expose COGEntry struct
    py::class_<COGEntry>(m, "COGEntry")
        .def(py::init<>())
        .def_readwrite("protein", &COGEntry::protein)
        .def_readwrite("organism", &COGEntry::organism)
        .def_readwrite("protein_name", &COGEntry::protein_name)
        .def_readwrite("length", &COGEntry::length)
        .def_readwrite("start", &COGEntry::start)
        .def_readwrite("end", &COGEntry::end)
        .def_readwrite("cog_id", &COGEntry::cog_id);
    
    // Expose AllEdges struct
    py::class_<AllEdges>(m, "AllEdges")
        .def(py::init<>())
        .def_readwrite("from_id", &AllEdges::from_id)
        .def_readwrite("to_id", &AllEdges::to_id);

    // Expose COGEdges struct
    py::class_<COGEdges>(m, "COGEdges")
        .def(py::init<>())
        .def_readwrite("cog_id", &COGEdges::cog_id)
        .def_readwrite("from_id", &COGEdges::from_id)
        .def_readwrite("to_id", &COGEdges::to_id);
    
    // Expose COGResults struct
    py::class_<COGResults>(m, "COGResults")
        .def(py::init<>())
        .def_readwrite("cog_entries", &COGResults::cog_entries)
        .def_readwrite("all_edges", &COGResults::all_edges)
        .def_readwrite("cog_edges", &COGResults::cog_edges);
    
    // Expose HitSet
    py::class_<HitSet>(m, "HitSet")
    .def(
        py::init([](double expect,
                    double overlap,
                    py::dict blast_dict,
                    py::dict prot2org_dict,
                    py::dict prot2len_dict,
                    py::dict all_lse_dict,
                    py::dict prot2lse_dict,
                    std::string cog_name,
                    int start_num,
                    py::dict query2subject_dict=py::none()) {
            // Convert blast hits
            auto blast_hits = convert_blast_hits(blast_dict);

            // Convert prot2org
            std::map<int, std::string> prot2org;
            for (const auto& item : prot2org_dict) {
                prot2org[item.first.cast<int>()] =
                    item.second.cast<std::string>();
            }

            // Convert prot2len
            std::map<int, int> prot2len;
            for (const auto& item : prot2len_dict) {
                prot2len[item.first.cast<int>()] =
                    item.second.cast<int>();
            }

            // Convert all_lse
            std::map<int, Lse> all_lse;
            for (const auto& item : all_lse_dict) {
                int lse_id = item.first.cast<int>();
                py::list prots = item.second.cast<py::list>();
                Lse lse_obj;
                for (auto p : prots) {
                    lse_obj.lse.insert(p.cast<int>());
                }
                all_lse[lse_id] = lse_obj;
            }

            // Convert prot2lse
            std::map<int, int> prot2lse;
            for (const auto& item : prot2lse_dict) {
                prot2lse[item.first.cast<int>()] =
                    item.second.cast<int>();
            }

            // Convert query2subject
            auto query2subject = convert_query2subject(query2subject_dict);

            auto cog_logger = spdlog::get("cogmaker");

            // Construct and return HitSet
            return HitSet(
                expect, overlap,
                blast_hits,
                prot2org, prot2len,
                all_lse, prot2lse,
                cog_name, start_num,
                query2subject,
                cog_logger
            );
        }),
        py::arg("expect"),
        py::arg("overlap"),
        py::arg("blast"),
        py::arg("prot2org"),
        py::arg("prot2len"),
        py::arg("all_lse"),
        py::arg("prot2lse"),
        py::arg("cog_name"),
        py::arg("start_num"),
        py::arg("query2subject")
    )
    // Enables calling of HitSet.insert() from Python
    .def("insert", [](HitSet& self, py::set lse_set) {
        std::set<int> cpp_set;
        for (auto item : lse_set) {
            cpp_set.insert(item.cast<int>());
        }
        self.insert(cpp_set);
    })
    // Enables calling of HitSet.makeCOGs() from Python
    .def("makeCOGs", [](HitSet& self) {
        return self.makeCOGs();
    });
};