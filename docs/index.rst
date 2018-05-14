#############################
Interstellar Dust Attenuation
#############################

``dust_attenuation`` is a python package to provide interstellar dust
attenuation curves.

Attenuation describe the effects of dust on a group of stars.  The effects
include in attenuation are dust absorption, dust scattering out of the
line-of-sight, and dust scattering into the line-of-sight.  In general,
attenuation models are used to model or correct the effects of dust on
observations of region of galaxies or global measurements of galaxies.

In contrast, dust extinction refers to the effects of dust on measurements
of a single star.  In this case, dust absorbs and scatters photons out of
the line-of-sight, but the contribution from scattered photons into the
line-of-sight is negligible.  For extinction models, see
the `dust_attenuation <https://github.com/karllark/dust_attenuation>`_ package.


This package is developed in the
`astropy affiliated package <http://www.astropy.org/affiliated/>`_
template and uses the
`astropy.modeling <http://docs.astropy.org/en/stable/modeling/>`_
framework.

User Documentation
==================

.. toctree::
   :maxdepth: 2

   Flavors of Models <dust_attenuation/model_flavors.rst>
   Fitting attenuation curves <dust_attenuation/fit_attenuation.rst>
   How to choose a model <dust_attenuation/choose_model.rst>
   References <dust_attenuation/references.rst>

Installation
============

.. toctree::
  :maxdepth: 2

  How to install <dust_attenuation/install.rst>

Quick Start
===========

Need to add.

Reporting Issues
================

If you have found a bug in ``dust_attenuation`` please report it by creating a
new issue on the ``dust_attenuation`` `GitHub issue tracker
<https://github.com/karllark/dust_attenuation/issues>`_.

Please include an example that demonstrates the issue sufficiently so that
the developers can reproduce and fix the problem. You may also be asked to
provide information about your operating system and a full Python
stack trace.  The developers will walk you through obtaining a stack
trace if it is necessary.

Contributing
============

Like the `Astropy`_ project, ``dust_attenuation`` is made both by and for its
users.  We accept contributions at all levels, spanning the gamut from
fixing a typo in the documentation to developing a major new feature.
We welcome contributors who will abide by the `Python Software
Foundation Code of Conduct
<https://www.python.org/psf/codeofconduct/>`_.

``dust_attenuation`` follows the same workflow and coding guidelines as
`Astropy`_.  The following pages will help you get started with
contributing fixes, code, or documentation (no git or GitHub
experience necessary):

* `How to make a code contribution <http://astropy.readthedocs.io/en/stable/development/workflow/development_workflow.html>`_

* `Coding Guidelines <http://docs.astropy.io/en/latest/development/codeguide.html>`_

* `Try the development version <http://astropy.readthedocs.io/en/stable/development/workflow/get_devel_version.html>`_

* `Developer Documentation <http://docs.astropy.org/en/latest/#developer-documentation>`_


For the complete list of contributors please see the `dust_extinction
contributors page on Github
<https://github.com/karllark/dust_attenuation/graphs/contributors>`_.

Reference API
=============

.. automodapi:: dust_attenuation.C00

.. automodapi:: dust_attenuation.WG00
