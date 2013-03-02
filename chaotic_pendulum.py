#!/usr/bin/env python

"""
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
"""

import numpy
import scipy.integrate
import matplotlib.pyplot, matplotlib.animation
import csv

from numpy import sin, cos, pi, array

"""
    Double pendulum constants
"""

M1 = 1
M2 = 1
L1 = 1
L2 = 1
G = 9.80665 # Gravitation, in m/s^2

"""
    Initial condition
"""

THETA1 = 90 # Initial value of theta1, in degrees
THETA2 = 90 # Initial value of theta1, in degrees
W1 = 2      # Initial angular velocity for first pendulum, in rad/s
W2 = 2      # Initial angular velocity for second pendulum, in rad/s

"""
    Integration constants
"""

T_MAX = 20  # End time of simulation, in seconds
DT = 0.05   # Integration step width, in seconds

def derivs(state, t):
    """
        State: 4-Tuple with following elements: (theta1, omega1, theta2, omega2)
        t:

        The equations can be found on http://www.physics.usyd.edu.au/~wheat/dpend_html/
        An image of the system can be found on http://bit.ly/YK88T5
    """

    t1 = state[0]   # Theta1, the angle of the first pendulum
    w1 = state[1]   # Omega1, the angular velocity for the first pendulum
    t2 = state[2]   # Theta2, the angle of the second pendulum
    w2 = state[3]   # Omega2, the angular velocity for the second pendulum

    # The solving

    delta = t2 - t1

    theta1_dot = w1
    theta2_dot = w2

    omega1_dot = ( M2 * L1 * w1 * w1 * sin(delta) * cos(delta) \
               + M2 *G * sin(t2) * cos(delta) \
               + M2 * L2 * w2 * w2 * sin(delta) \
               - (M1 + M2 ) * G *sin(t1) ) \
               / ( (M1 + M2 ) * L1 - M2 * L1 * cos(delta) * cos(delta) )

    omega2_dot = ( - M2 * L2 * w2 * w2  * sin(delta) * cos(delta) \
               + ( M1 + M2 ) \
               * (G * sin(t1) * cos(delta) - L1 * w1 * w1 * sin(delta) - G * sin(t2))) \
               / ( (M1 + M2 ) * L2 - M2 * L2 * cos(delta) * cos(delta) )

    return (theta1_dot, omega1_dot, theta2_dot, omega2_dot)

def get_positions(y):
    """ Calculates the x- and y-coordinates of the pendulum from the related angles
    """

    x1 = L1*sin(y[:,0])
    y1 = -L1*cos(y[:,0])

    x2 = L2*sin(y[:,2]) + x1
    y2 = -L2*cos(y[:,2]) + y1

    return (x1, y1, x2, y2)

def dump_csv(t, y, filename):
    """ Dumps the time, positions, angles and velocitys in a csv file
        t: Time data of the system
        y: state array of the system+
        filename: name of the file to save to
    """

    with open(filename, 'wb') as f:
        fieldnames = 't theta1 omega1 theta2 omega2 x1 y1 x2 y2'.split()
        writer = csv.DictWriter(f, fieldnames, delimiter=';')
        writer.writeheader()

        x1, y1, x2, y2 = get_positions(y)

        for i in range( len(t) ):
            dt = t[i]
            t1, w1, t2, w2 = y[i]

            d = {
                "t" : dt,
                "theta1" : t1,
                "omega1" : w1,
                "theta2" : t2,
                "omega2" : w2,
                "x1" : x1[i],
                "y1" : y1[i],
                "x2" : x2[i],
                "y2" : y2[i]
            }
            writer.writerow(d)

def movietime(t, y):
    x1, y1, x2, y2 = get_positions(y)
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(-2, 2), ylim=(-2, 2))
    ax.grid()
    line, = ax.plot([], [], 'o-', lw=2)
    time_template = 'time = %.1fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

    def init():
        line.set_data([], [])
        time_text.set_text('')
        return line, time_text

    def animate(i):
        thisx = [0, x1[i], x2[i]]
        thisy = [0, y1[i], y2[i]]

        line.set_data(thisx, thisy)
        time_text.set_text(time_template%(i*DT))
        return line, time_text


    ani = matplotlib.animation.FuncAnimation(fig, animate, numpy.arange(1, len(y)), interval=25, blit=True, init_func=init)
    Writer = matplotlib.animation.writers['ffmpeg']
    writer = Writer(fps=15)
    ani.save('double_pendulum.mp4', writer=writer)

if __name__ == '__main__':
    state = numpy.array([THETA1, W1, THETA2, W2])*pi/180.
    t = numpy.arange(0.0, T_MAX, DT)
    y = scipy.integrate.odeint(derivs, state, t)

    # dump_csv(t, y, 'pendulum.csv')
    movietime(t,y)








