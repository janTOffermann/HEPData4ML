/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "display/Utils.h"
#include <vector>
#include "TColor.h"

using namespace std;

namespace Display{

vector<EColor> GetRainbowColorsCustom(Int_t nColors,
                                          Float_t saturation,
                                          Float_t brightness) {
    vector<EColor> colors;

    if (nColors <= 0) return colors;

    colors.reserve(nColors);

    for (Int_t i = 0; i < nColors; ++i) {
        Float_t hue = (300.0 * i) / (nColors - 1);

        Float_t r, g, b;
        TColor::HSV2RGB(hue, saturation, brightness, r, g, b);

        Int_t colorIndex = TColor::GetFreeColorIndex();
        TColor* newColor = new TColor(colorIndex, r, g, b);

        colors.push_back(static_cast<EColor>(colorIndex));
    }
    return colors;
}

}
//------------------------------------------------------------------------------
