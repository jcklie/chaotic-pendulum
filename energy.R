#!/bin/env R

"
    Simulation of a double pendulum. For more information, feel free
    to visit the home: https://github.com/Rentier/chaotic-pendulum

    Copyright (C) 2013 Jan-Christoph Klie

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"

data <- read.table("double_pendulum.csv", header=T,sep=";")

plot(x=data[ , "t"], y=data[ , "energy"], xlab="t [s]", ylab="E [J]")
