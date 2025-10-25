""" fd3d_5_2.py: Interactive

Use interactive widgets to view output data from fdtd_5_1
"""
from collections import namedtuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import RadioButtons

plt.rcParams['toolbar'] = 'None'


def main():
    plot_parameters = load_plot_parameters()

    figure, current_ax, fdtd_plot, bessel_plot = plot_fig_4_7(
        plot_parameters[0])
    plt.subplots_adjust(left=0.33)

    # Create an instance of the class 'Controller'
    radio_controller = Controller(figure, current_ax, plot_parameters,
                                  fdtd_plot, bessel_plot)
    plt.show()


class Controller:
    """ This class creates the controller for the interactive widgets.

    Attributes:
        figure: figure
        current_ax: axes for the figure
        plot_parameters: data saved from fd3d_5_1 (frequency, fdtd
                    amplitude, Bessel amplitude, and x-axis data)
        fdtd_plot = fdtd data plot
        bessel_plot = Bessel data plot
    """

    def __init__(self, figure, current_ax, plot_parameters, fdtd_plot,
                 bessel_plot):
        self.figure = figure
        self.current_ax = current_ax
        self.plot_parameters = plot_parameters
        self.fdtd_plot = fdtd_plot
        self.bessel_plot = bessel_plot

        self.selected_freq_label = \
            plot_parameters[0].frequency_label  # Initial frequency

        # This sets up the radio buttons
        radio_axes = plt.axes([0.03, 0.55, 0.15, 0.20])  # Radio button axis
        self.radio = RadioButtons(
            ax=radio_axes,
            labels=[output.frequency_label for output in plot_parameters]
        )
        self.radio.on_clicked(self.on_radio_select)

    def on_radio_select(self, label_selected):
        """On radio select, the new frequency is selected
        and plots are redrawn"""
        self.selected_freq_label = label_selected
        self.redraw()

    def redraw(self):
        """Redraw the figure based on current_freq"""
        for plot_parameter in self.plot_parameters:
            if self.selected_freq_label == plot_parameter.frequency_label:
                self.fdtd_plot.set_data(plot_parameter.fdtd_location,
                                        plot_parameter.fdtd_amplitude)
                self.bessel_plot.set_data(plot_parameter.bessel_location,
                                          plot_parameter.bessel_amplitude)
        self.current_ax.relim()
        self.current_ax.autoscale_view(True, True, True)
        plt.savefig('figure_5_2')
        plt.draw()


PlotParameters = namedtuple('PlotParameters',
                            ['frequency', 'frequency_label',
                             'fdtd_location', 'fdtd_amplitude',
                             'bessel_location', 'bessel_amplitude'])


def load_plot_parameters():
    """ Load the arrays containing the outputs from fd3d_5_1 """

    fdtd_x_axis = np.load('fdtd_x_axis.npy')
    fdtd_amp = np.load('fdtd_amp.npy')
    bessel_x_axis = np.load('bessel_x_axis.npy')
    bessel_amp = np.load('bessel_amp.npy')
    frequencies = np.load('frequencies.npy')

    return [
        PlotParameters(
            frequency=frequency,
            frequency_label="{} MHz".format(int(frequency / 1e6)),
            fdtd_location=fdtd_x_axis,
            fdtd_amplitude=fdtd_amp[i, :],
            bessel_location=bessel_x_axis,
            bessel_amplitude=bessel_amp[i, :]
        )
        for i, frequency in enumerate(frequencies)
    ]


def plot_fig_4_7(plot_parameter):
    """Plot Fig. 5.2
    This is used to view one frequency at a time"""

    plt.rcParams['font.size'] = 12
    plt.rcParams['grid.color'] = 'gray'
    plt.rcParams['grid.linestyle'] = 'dotted'
    figure, current_ax = plt.subplots(figsize=(8, 3))

    fdtd_plot, = current_ax.plot(plot_parameter.fdtd_location,
                                 plot_parameter.fdtd_amplitude,
                                 color='k',
                                 linewidth=1)
    bessel_plot, = current_ax.plot(plot_parameter.bessel_location,
                                   plot_parameter.bessel_amplitude,
                                   'ko',
                                   mfc='none',
                                   linewidth=1)
    plt.xlabel('cm')
    plt.ylabel('Amplitude')
    plt.xticks(np.arange(-5, 10, step=5))
    plt.xlim(-9, 9)

    plt.tight_layout()
    return figure, current_ax, fdtd_plot, bessel_plot


if __name__ == '__main__':
    main()
