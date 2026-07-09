/* ----------------------------------------------------------------------
* I-SIMPA (https://i-simpa.univ-gustave-eiffel.fr). This file is part of I-SIMPA.
*
* I-SIMPA is a GUI for 3D numerical sound propagation modelling dedicated
* to scientific acoustic simulations.
* Copyright (C) UMRAE, CEREMA, Univ Gustave Eiffel - Judicael Picaut, Nicolas Fortin
*
* I-SIMPA is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or
* (at your option) any later version.
*
* I-SIMPA is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software Foundation,
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA or
* see <http://ww.gnu.org/licenses/>
*
* For more information, please consult: <https://i-simpa.univ-gustave-eiffel.fr> or
* send an email to contact@noise-planet.org
*
* ----------------------------------------------------------------------*/

#ifndef __AMBISONICS_EXPORT__
#define __AMBISONICS_EXPORT__

#include <sppsTypes.h>
#include <vector>

typedef dvec3 veci_t;

/**
 * @file ambisonics_export.h
 * @brief First-order B-format (WXYZ, FuMa) impulse response export for punctual receivers
 */
namespace AmbisonicsExport
{
	/**
	 * One frequency band's per-time-bin energy and energy-weighted direction, as tracked
	 * by the calculation core (either combined across all sources, or for a single source).
	 */
	struct t_BandContribution
	{
		double centerFreq;					/*!< Band center frequency (Hz) */
		const l_decimal* energyPerBin;		/*!< Energy received in this band, per time step (size nbTimeStep) */
		const veci_t* directionPerBin;		/*!< Energy-weighted direction vector, per time step (size nbTimeStep) */
	};

	/**
	 * Synthesizes and writes a 4-channel (W,X,Y,Z) 32-bit float B-format impulse response WAV file
	 * from per-band energy/direction contributions tracked by the calculation core.
	 */
	bool WriteBFormatImpulseResponse(const CoreString& wavFilePath, const std::vector<t_BandContribution>& bands,
		int nbTimeStep, double timeStep, unsigned int randomSeed);
}

#endif // __AMBISONICS_EXPORT__
