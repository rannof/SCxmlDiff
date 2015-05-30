# SCxmlDiff
###Compare two SeisComp3 event parameter xml files.
Created by Ran Novitsky Nof (ran.rnof.info), 2015 @ BSL.

#### DEPENDENCIES:
-   numpy - http://numpy.org
-   matplotlib - http://matplotlib.org
-   PyQt4 - http://www.riverbankcomputing.com
-   pyproj - http://jswhit.github.io/pyproj
-   SeisComP3 - http://www.seiscomp3.org

#### USAGE:
<pre>
SCxmlDiff.py [-h] [-i file [file ...]] [-r file [file ...]]
                    [-b bounds bounds bounds bounds]

optional arguments:
  -h, --help            show this help message and exit
  -i file [file ...]    input E2 result files (xml or log)
  -r file [file ...]    input reference files (xml or csv) xml - Seiscomp3 xml
                        file csv - comma separated value (GII DB forat
                        [epiid,ml,typ,lat,lon,abs_ot_d,abs_ot_t].
  -b bounds bounds bounds bounds
                        Region bounding box (west east south north)
</pre>
#### LICENSE:
  SCxmlDiff is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
