
#include "plotting.hpp"

namespace plt = matplotlibcpp;

void plotting::plot_evolution(const std::vector<integrator::step>& results, const std::vector<integrator::Event>& events, std::vector<int> plot_components, std::vector<std::string> labels) {
    assert(results.size() > 0);
    assert(plot_components.size() == labels.size());
    std::vector<double> r;
    r.reserve(results.size());
    for(auto it = results.begin(); it != results.end(); ++it)
        r.push_back(it->first);

    for(int i = 0; i < plot_components.size(); i++) {
        std::vector<double> y;
        y.reserve(results.size());
        int index = plot_components[i];
        assert(index >= 0 && index < results[0].second.size());
        for(auto it = results.begin(); it != results.end(); ++it)
            y.push_back(it->second[index]);
        plt::plot(r, y, {{"label", labels[i]}});
    }
}

