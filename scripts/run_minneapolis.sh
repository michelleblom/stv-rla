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


python3.9 audit_general.py -d data/Minneapolis/MPLS-2009-BET_2Seat_ParsedMB.txt -s 2  -outcome data/Minneapolis/MPLS-2009-BET.otc -quota 10696 -voters 47857 -log data/Minneapolis/MPLS-2009-BET.log

python3.9 audit_general.py -d data/Minneapolis/MPLS-2013-BET_2Seat_ParsedMB.txt -s 2  -outcome data/Minneapolis/MPLS-2013-BET.otc -quota 16286 -voters 48855 -log data/Minneapolis/MPLS-2013-BET.log

python3.9 audit_general.py -d data/Minneapolis/MPLS-2017-BET_2Seat_ParsedMB.txt -s 2  -outcome data/Minneapolis/MPLS-2017-BET.otc -quota 23232 -voters 69694 -log data/Minneapolis/MPLS-2017-BET.log -gelim 2

python3.9 audit_general.py -d data/Minneapolis/MPLS-2021-BET_2Seat_ParsedMB.txt -s 2  -outcome data/Minneapolis/MPLS-2021-BET.otc -quota 31876 -voters 145337 -log data/Minneapolis/MPLS-2021-BET.log -gelim 2



