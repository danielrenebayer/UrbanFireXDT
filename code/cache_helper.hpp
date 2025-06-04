/**
 * cache_helper.hpp
 *
 * Header-only file containing the stuff for caching allocated profile IDs
 **/

#pragma once

// This header files defines the following classes
class HPProfileIDCache;
class PVProfileIDCache;

#include <boost/json.hpp>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <optional>
#include <string>
#include <mutex>

namespace json = boost::json;

/**
 * This is the singelton class for storing and querying the cache file for
 * the heat pump profiles allocated to a given building referenced by the unitID.
 */
class HPProfileIDCache {
    private:
        HPProfileIDCache() = default;
        ~HPProfileIDCache() = default;

    public:
        HPProfileIDCache(HPProfileIDCache const&) = delete;
        void operator=(HPProfileIDCache const&)   = delete;

        /**
         * Static method for getting the only instance of the singleton class
         */
        static HPProfileIDCache& GetInstance() {
            static HPProfileIDCache instance_;  // thread-safe since C++11, instantiated on first use
            return instance_;
        }

        /**
         *  Set the filename for the cache file
         **/
        void setCacheFilename(const std::string& filename) {
            std::lock_guard<std::mutex> lock(mutex_);
            cache_filename_ = filename;
            load_from_file();
        }

        /**
         * Update or insert an entry
         **/
        void updateCache(size_t unitID, size_t profileID) {
            std::lock_guard<std::mutex> lock(mutex_);
            cache_updated = true;
            cache_[unitID] = profileID;
            //save_to_file();
        }

        /**
         * Read an entry. Returns std::nullopt if not found
         **/
        std::optional<size_t> readCache(size_t unitID) {
            auto it = cache_.find(unitID);
            if (it != cache_.end()) {
                return it->second;
            }
            return std::nullopt;
        }

        /**
         * Saves the current content of the cache file.
         */
        void saveCacheFile() {
            std::lock_guard<std::mutex> lock(mutex_);
            save_to_file();
        }

    private:
        std::string cache_filename_;
        std::unordered_map<size_t, size_t> cache_;
        bool cache_updated = false;
        std::mutex mutex_;

        void save_to_file() {
            if (!cache_updated)
                return;

            if (cache_filename_.empty()) {
                throw std::runtime_error("HP cache filename is not set. Call set_cache_filename() before saving.");
            }

            json::array jarray;
            for (const auto& [unitID, profileID] : cache_) {
                jarray.emplace_back(json::array{unitID, profileID});
            }

            std::ofstream ofs(cache_filename_);
            if (ofs) {
                ofs << json::serialize(jarray);
            }

            std::cout << "Heat pump profile allocation written to cache file '" << cache_filename_ << "'\n";
        }

        void load_from_file() {
            cache_.clear();
            std::ifstream ifs(cache_filename_);
            if (!ifs) return;

            std::stringstream buffer;
            buffer << ifs.rdbuf();

            boost::system::error_code ec;
            json::value jv = json::parse(buffer.str(), ec);
            if (ec) {
                throw std::runtime_error("Failed to parse JSON cache file: " + ec.message());
            }
            const json::array& jarray = jv.as_array();

            for (const auto& item : jarray) {
                const json::array& pair = item.as_array();
                size_t unitID;
                size_t profileID;
                bool skip = false;

                if (pair.size() != 2) {
                    skip = true;
                }
                if (!skip && pair[0].is_int64()) {
                    unitID = static_cast<size_t>(pair[0].as_int64());
                } else if (!skip && pair[0].is_uint64()) {
                    unitID = static_cast<size_t>(pair[0].as_uint64());
                } else {
                    skip = true;
                }
                if (!skip && pair[1].is_int64()) {
                    profileID = static_cast<size_t>(pair[1].as_int64());
                } else if (!skip && pair[1].is_uint64()) {
                    profileID = static_cast<size_t>(pair[1].as_uint64());
                } else {
                    skip = true;
                }

                if (skip) {
                    std::cerr << "Skipping malformed cache entry: " << json::serialize(pair) << "\n";
                    continue;
                }

                cache_[unitID] = profileID;
            }

            std::cout << "Heat pump profile allocation read from cache file '" << cache_filename_ << "'\n";
        }
};



/**
 * This is the singleton class for storing and querying the cache file for
 * the PV profiles allocated to a given location referenced by the locationID and roof orientation.
 */
class PVProfileIDCache {
    private:
        PVProfileIDCache() = default;
        ~PVProfileIDCache() = default;

    public:
        PVProfileIDCache(PVProfileIDCache const&) = delete;
        void operator=(PVProfileIDCache const&)   = delete;

        /**
         * Static method for getting the only instance of the singleton class
         */
        static PVProfileIDCache& GetInstance() {
            static PVProfileIDCache instance_;  // thread-safe since C++11, instantiated on first use
            return instance_;
        }

        /**
         *  Set the filename for the cache file
         **/
        void setCacheFilename(const std::string& filename) {
            std::lock_guard<std::mutex> lock(mutex_);
            cache_filename_ = filename;
            load_from_file();
        }

        /**
         * Update or insert an entry
         **/
        void updateCache(size_t locationID, const std::string& roof_orientation, size_t profileID) {
            std::lock_guard<std::mutex> lock(mutex_);
            cache_updated = true;
            cache_[locationID][roof_orientation] = profileID;
            //save_to_file();
        }

        /**
         * Read an entry. Returns std::nullopt if not found
         **/
        std::optional<size_t> readCache(size_t locationID, const std::string& roof_orientation) {
            std::lock_guard<std::mutex> lock(mutex_);
            auto loc_it = cache_.find(locationID);
            if (loc_it != cache_.end()) {
                auto orient_it = loc_it->second.find(roof_orientation);
                if (orient_it != loc_it->second.end()) {
                    return orient_it->second;
                }
            }
            return std::nullopt;
        }

        /**
         * Saves the current content of the cache file.
         */
        void saveCacheFile() {
            std::lock_guard<std::mutex> lock(mutex_);
            save_to_file();
        }

    private:
        std::string cache_filename_;
        std::unordered_map<size_t, std::unordered_map<std::string, size_t>> cache_; // locationID -> map< orientation string, profileID >
        bool cache_updated = false;
        std::mutex mutex_;

        void save_to_file() {
            if (!cache_updated)
                return;

            if (cache_filename_.empty()) {
                throw std::runtime_error("PV cache filename is not set. Call set_cache_filename() before saving.");
            }

            json::array jarray;
            for (const auto& [locationID, orientation_map] : cache_) {
                for (const auto& [orientation, profileID] : orientation_map) {
                    jarray.emplace_back(json::array{locationID, orientation, profileID});
                }
            }

            std::ofstream ofs(cache_filename_);
            if (ofs) {
                ofs << json::serialize(jarray);
            }

            std::cout << "PV profile allocation written to cache file '" << cache_filename_ << "'\n";
        }

        void load_from_file() {
            cache_.clear();
            std::ifstream ifs(cache_filename_);
            if (!ifs) return;

            std::stringstream buffer;
            buffer << ifs.rdbuf();

            boost::system::error_code ec;
            json::value jv = json::parse(buffer.str(), ec);
            if (ec) {
                throw std::runtime_error("Failed to parse JSON cache file: " + ec.message());
            }

            const json::array& jarray = jv.as_array();
            for (const auto& item : jarray) {
                const json::array& entry = item.as_array();
                size_t locationID;
                std::string roof_orientation;
                size_t profileID;
                bool skip = false;

                if (entry.size() != 3) {
                    skip = true;
                }
                if (!skip && entry[0].is_int64()) {
                    locationID = static_cast<size_t>(entry[0].as_int64());
                } else if (!skip && entry[0].is_uint64()) {
                    locationID = static_cast<size_t>(entry[0].as_uint64());
                } else {
                    skip = true;
                }
                if (!skip && entry[1].is_string()) {
                    roof_orientation = entry[1].as_string().c_str();
                } else {
                    skip = true;
                }
                if (!skip && entry[2].is_int64()) {
                    profileID = static_cast<size_t>(entry[2].as_int64());
                } else if (!skip && entry[2].is_uint64()) {
                    profileID = static_cast<size_t>(entry[2].as_uint64());
                } else {
                    skip = true;
                }

                if (skip) {
                    std::cerr << "Skipping malformed cache entry: " << json::serialize(entry) << "\n";
                    continue;
                }

                cache_[locationID][roof_orientation] = profileID;
            }

            std::cout << "PV profile allocation read from cache file '" << cache_filename_ << "'\n";
        }
};


