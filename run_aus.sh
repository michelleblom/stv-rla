#
#    Copyright (C) 2021  Michelle Blom
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

python3 audit_general.py -d data/FedAus/federal_2016_ACT.stv -outcome data/FedAus/federal_2016_ACT.otc -quota 84923 -voters 260521 -e 0.002 -r 0.10 -log "data/FedAus/federal_2016_ACT.0.10.0.002.log" 

python3 audit_general.py -d data/FedAus/federal_2019_ACT.stv -outcome data/FedAus/federal_2019_ACT.otc -quota 90078 -voters 276651 -e 0.002 -r 0.10 -log "data/FedAus/federal_2019_ACT.0.10.0.002.log" 

python3 audit_general.py -d data/FedAus/federal_2016_NT.stv -outcome data/FedAus/federal_2016_NT.otc -quota 34010 -voters 105539 -e 0.002 -r 0.10 -log "data/FedAus/federal_2016_NT.0.10.0.002GEN.log" 

python3 audit_general.py -d data/FedAus/federal_2019_NT.stv -outcome data/FedAus/federal_2019_NT.otc -quota 35010 -voters 108994 -e 0.002 -r 0.10 -log "data/FedAus/federal_2019_NT.0.10.0.002GEN.log" 


