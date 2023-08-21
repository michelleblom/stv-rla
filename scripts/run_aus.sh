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

python3.9 audit_general.py -d data/FedAus/federal_2016_ACT.stv -outcome data/FedAus/federal_2016_ACT.otc -quota 84923 -voters 260521 -e1 0.002 -e2 0 -r 0.10 -log "data/FedAus/federal_2016_ACT.0.10.0.002.log" 

python3.9 audit_general.py -d data/FedAus/federal_2019_ACT.stv -outcome data/FedAus/federal_2019_ACT.otc -quota 90078 -voters 276651 -e1 0.002 -e2 0 -r 0.10 -log "data/FedAus/federal_2019_ACT.0.10.0.002.log" 

python3.9 audit_general.py -d data/FedAus22/2022ACT.json -outcome data/FedAus22/2022ACT.otc_D2 -quota 90078 -voters 276651 -e1 0.002 -e2 0 -r 0.10 -log "data/FedAus/federal_2022_ACT.0.10.0.002.log" 



python3.9 audit_general.py -d data/FedAus/federal_2016_NT.stv -outcome data/FedAus/federal_2016_NT.otc -quota 34010 -voters 102027 -e1 0.002 -e2 0 -r 0.10 -log "data/FedAus/federal_2016_NT.0.10.0.002.log" 

python3.9 audit_general.py -d data/FedAus/federal_2019_NT.stv -outcome data/FedAus/federal_2019_NT.otc -quota 35010 -voters 105027 -e1 0.002 -e2 0 -r 0.10 -log "data/FedAus/federal_2019_NT.0.10.0.002.log" 

