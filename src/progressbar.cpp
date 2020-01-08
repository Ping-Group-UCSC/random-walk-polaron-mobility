#include "progressbar.h"

bool ProgressBar::is_disable = true;

ProgressBar::ProgressBar(std::string const& name, size_t max, size_t width):
    name(name), max(max), width(width), last(max*2)
{
    if (ProgressBar::is_disable) return;
    std::cout << std::endl;
}

void ProgressBar::Update(size_t val) {
    if (ProgressBar::is_disable) return;
    std::cout << "\r" << name << "[";
    int percent = 100 * val / max;
    //Do not update if percent does not change
    if (100 * val / max != 100 * last / max) {
        int i = 0;
        for (; i < width * val / max ; i++){
            std::cout << "=";
        }
        for (; i < width ; i++){
            std::cout << " ";
        }
        std::cout << "] " <<  val << "/" << max << " " << 100 * val / max << "%"  ;
    }
    last = val;
    std::cout.flush();
}

void ProgressBar::Add(size_t n){
    if (ProgressBar::is_disable) return;
    Update(last + n);
}

void ProgressBar::Clear() {
    if (ProgressBar::is_disable) return;
    std::cout << "\x1b[2K\r\x1b[A" ;
}

