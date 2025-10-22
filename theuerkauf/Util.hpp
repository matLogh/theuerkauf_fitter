#pragma once

#ifndef __Util_h__
#define __Util_h__

#include <mutex>
#include <sstream>
#include <string>

#include <TAxis.h>
#include <TColor.h>
#include <TH1.h>

inline std::string GetFuncUniqueName(const char *prefix, void *ptr)
{
    static int ___num = 0;
    static std::mutex __getfunctuniquename__mtx;
    std::lock_guard<std::mutex> lock(__getfunctuniquename__mtx);

    std::ostringstream name;
    name << prefix << "_" << ptr << "_" << ++___num;

    return name.str();
}

inline double TH1IntegrateWithPartialBins(const TH1 *spec, const double xmin, const double xmax)
{
    const TAxis *axis = spec->GetXaxis();
    const int bmin = axis->FindBin(xmin);
    const int bmax = axis->FindBin(xmax);
    double integral = spec->Integral(bmin, bmax);
    integral -= spec->GetBinContent(bmin) * (xmin - axis->GetBinLowEdge(bmin)) / axis->GetBinWidth(bmin);
    integral -= spec->GetBinContent(bmax) * (axis->GetBinUpEdge(bmax) - xmax) / axis->GetBinWidth(bmax);
    return integral;
}

/// @brief check if the hex color code is in a valid format, such as "#c0c0c0"
/// @param input
/// @return
inline bool isValidHexColor(const std::string &input)
{
    if (input.size() != 7 || input[0] != '#')
    {
        return false; // Hex color codes are of the form "#RRGGBB"
    }

    // Convert any uppercase letters to lowercase
    std::string lowerInput = input;
    std::transform(lowerInput.begin(), lowerInput.end(), lowerInput.begin(), ::tolower);

    for (size_t i = 1; i < lowerInput.size(); ++i)
    {
        char c = lowerInput[i];
        if (!isxdigit(c) || (c >= 'g' && c <= 'z'))
        {
            return false; // Must be a valid hexadecimal digit (0-9, a-f)
        }
    }

    return true;
}

/// @brief Thread-safe way to get a color from a predifned list. It is also possible to permanently add new color to the
/// list using the input parameter
/// @param hex hex color code, must be in format "#c0c0c0"
/// @param reset reset color index from 0
/// @return color_t used by ROOT (but it is just Int_t)
inline Color_t GetColor(int get_color_index = -1, std::string hex_color = "")
{
    static std::mutex __getcolor__mtx;
    std::lock_guard<std::mutex> lock(__getcolor__mtx);
    static unsigned int ___color_index = 0;
    static std::vector<std::string> ___color_array{
        "#2196f3", "#f44336", "#3f51b5", "#4caf50", "#ff9800", "#000000", "#e91e63", "#8bc34a", "#ff5722", "#795548",
        "#607d8b", "#9c27b0", "#00bcd4", "#ffeb3b", "#673ab7", "#03a9f4", "#009688", "#cddc39", "#ffc107", "#9e9e9e"};
    if (!hex_color.empty() && isValidHexColor(hex_color))
        ___color_array.emplace_back(hex_color);
    if (get_color_index > 0 && get_color_index < static_cast<int>(___color_array.size()))
        ___color_index = get_color_index;

    unsigned int index = ___color_index;
    ___color_index = ___color_index + 1 == ___color_array.size() ? 0 : ___color_index + 1;
    Color_t c = TColor::GetColor(___color_array[index].c_str());
    return c;
}

/// @brief Clear all functions attached to a histogram
/// @param h
inline void ClearListOfFunctions(TH1 *h)
{
    if (h == nullptr)
        return;
    TList *list = h->GetListOfFunctions();
    if (list != nullptr)
    {

        TIter iter(list);
        TObject *obj = nullptr;
        while ((obj = iter()))
        {
            list->Remove(obj);
            delete obj;
        }
    }
}

#endif // __Util_h__