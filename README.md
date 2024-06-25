# CPOLYMESH

## Table of Contents

- [Description](#description)
- [Installation](#installation)
- [Usage](#usage)
- [License](#license)
- [Authors](#authors)
- [Thank you to](#thanks_to)
- [References](#references)

## Description 

CPOLYMESH provide Python codes for polynomial approximation on complex compact sets with connected complement, by Chebyshev-like admissible polynomial meshes on boundaries with piecewise polynomial parametrization. Such meshes have lower cardinality with respect to those previously known. They are used for polynomial least-squares, for the extraction of extremal interpolation sets of Fekete and Leja type,
as well as for the computation of the uniform norms (Lebesgue constants) of polynomial
projection operators.

## Installation

### Using Anaconda
1. **Install Anaconda**:
   - Download and install Anaconda from the [official website](https://www.anaconda.com/products/individual).

2. **Clone the Repository**:
   ```bash
   git clone https://github.com/DimitriKenne/your_project.git
   cd your_project

3. Launch Spyder or VSCode:
  - If using Spyder:
    - Open Spyder from Anaconda Navigator.
  - If using VSCode:
    - Open VSCode and select the Anaconda Python interpreter.
    - Open the Command Palette (Ctrl+Shift+P) and type Python: Select Interpreter.
    - Select the interpreter that corresponds to the Anaconda distribution.

### Direct Installation Without Anaconda
If users prefer not to use Anaconda, then you can install the dependencies directly using pip:

1. Clone the Repository:
    git clone https://github.com/DimitriKenne/your_project.git
    cd your_project


2. Install Dependencies:
    pip install numpy scipy matplotlib



## Usage

### demo:
 In order to run the demo you need to open the file demo.py, scroll till its bottom, there provide the following information:
    - The degree of interpolation, deg, greater than 0

        Example: deg = 2

    - Choose the domain by specifying a number between 0 and 19 inside define_domaine. 

        Example: domain = define_domain(0)

 and finally run the whole file.


 The demo does:
    - compute admissible meshes (AM)

    - compute interpolations points: discrete Leja points (DLP), pseudo Leja points (PLP), approximate Fekete points (AFK)

    - approximate the Lebesgue constant for each of these points as well as for the dicrete least squares approximation (lsqp). The lsqp is performed on the whole AM of degree deg.

    - plot each interpolation points set + AM of degree deg on separate figures.

    - plot the Lebesgue constants on one figure.

After running the file demo.py for the first time, you can also type demo(2,define_domain(0)) directly from the terminal in order to run the demo with deg=2 and domain Num. 0. 

Example

>> demo(2,define_domain(0))

 Demo on computing Lebesgue constant

================================================================================
domain                                   circle of center 0 and radius 1         
Interp. nodes                            Discrete Leja points                    
degree                                   1                                       
Lebesgue Constant AM                     1.3870398453221475                      
Relative err. est.                       0.08239220029239402                     
Absolute err. est.                       0.1142812647493136                      
Mesh Leb. const. degree                  1                                       
Mesh Leb. const. factor                  4                                       
Mesh Leb. const. constant                1.082392200292394                       
Card. mesh Leb. const.                   8                                       
number of interpolation points           2                                       
Mesh interp. nodes degree                1                                       
Mesh interp. nodes factor                2                                       
Card mesh for extracting the nodes       4                                       
================================================================================
================================================================================
domain                                   circle of center 0 and radius 1         
Interp. nodes                            Pseudo Leja points                      
degree                                   1                                       
Lebesgue Constant AM                     1.3870398453221475                      
Relative err. est.                       0.08239220029239402                     
Absolute err. est.                       0.1142812647493136                      
Mesh Leb. const. degree                  1                                       
Mesh Leb. const. factor                  4                                       
Mesh Leb. const. constant                1.082392200292394                       
Card. mesh Leb. const.                   8                                       
number of interpolation points           2                                       
Mesh interp. nodes degree                1                                       
Mesh interp. nodes factor                2                                       
Card mesh for extracting the nodes       4                                       
================================================================================
================================================================================
domain                                   circle of center 0 and radius 1         
Interp. nodes                            Approximate Fekete points               
degree                                   1                                       
Lebesgue Constant AM                     1.3870398453221475                      
Relative err. est.                       0.08239220029239402                     
Absolute err. est.                       0.1142812647493136                      
Mesh Leb. const. degree                  1                                       
Mesh Leb. const. factor                  4                                       
Mesh Leb. const. constant                1.082392200292394                       
Card. mesh Leb. const.                   8                                       
number of interpolation points           2                                       
Mesh interp. nodes degree                1                                       
Mesh interp. nodes factor                2                                       
Card mesh for extracting the nodes       4                                       
================================================================================
================================================================================
domain                                   circle of center 0 and radius 1         
Interp. nodes                            Discrete least square approximation     
degree                                   1                                       
Lebesgue Constant AM                     1.281457723870753                       
Relative err. est.                       0.08239220029239402                     
Absolute err. est.                       0.10558212145139444                     
Mesh Leb. const. degree                  1                                       
Mesh Leb. const. factor                  4                                       
Mesh Leb. const. constant                1.082392200292394                       
Card. mesh Leb. const.                   8                                       
number of interpolation points           4                                       
Mesh interp. nodes degree                1                                       
Mesh interp. nodes factor                2                                       
Card mesh for extracting the nodes       4                                       
================================================================================

 Demo on computing Lebesgue constant

================================================================================
domain                                   circle of center 0 and radius 1         
Interp. nodes                            Discrete Leja points                    
degree                                   2                                       
Lebesgue Constant AM                     2.9615705608064604                      
Relative err. est.                       0.08239220029239402                     
Absolute err. est.                       0.2440103148260236                      
Mesh Leb. const. degree                  2                                       
Mesh Leb. const. factor                  4                                       
Mesh Leb. const. constant                1.082392200292394                       
Card. mesh Leb. const.                   16                                      
number of interpolation points           3                                       
Mesh interp. nodes degree                2                                       
Mesh interp. nodes factor                2                                       
Card mesh for extracting the nodes       8                                       
================================================================================
================================================================================
domain                                   circle of center 0 and radius 1         
Interp. nodes                            Pseudo Leja points                      
degree                                   2                                       
Lebesgue Constant AM                     3.1231886753400913                      
Relative err. est.                       0.08239220029239402                     
Absolute err. est.                       0.25732638688955756                     
Mesh Leb. const. degree                  2                                       
Mesh Leb. const. factor                  4                                       
Mesh Leb. const. constant                1.082392200292394                       
Card. mesh Leb. const.                   16                                      
number of interpolation points           3                                       
Mesh interp. nodes degree                2                                       
Mesh interp. nodes factor                2                                       
Card mesh for extracting the nodes       8                                       
================================================================================
================================================================================
domain                                   circle of center 0 and radius 1         
Interp. nodes                            approximate Fekete points               
degree                                   2                                       
Lebesgue Constant AM                     1.9155346724445566                      
Relative err. est.                       0.08239220029239402                     
Absolute err. est.                       0.1578251163990773                      
Mesh Leb. const. degree                  2                                       
Mesh Leb. const. factor                  4                                       
Mesh Leb. const. constant                1.082392200292394                       
Card. mesh Leb. const.                   16                                      
number of interpolation points           3                                       
Mesh interp. nodes degree                2                                       
Mesh interp. nodes factor                2                                       
Card mesh for extracting the nodes       8                                       
================================================================================
================================================================================
domain                                   circle of center 0 and radius 1         
Interp. nodes                            Discrete least square approximation     
degree                                   2                                       
Lebesgue Constant AM                     1.4339125628626894                      
Relative err. est.                       0.08239220029239402                     
Absolute err. est.                       0.11814321108116274                     
Mesh Leb. const. degree                  2                                       
Mesh Leb. const. factor                  4                                       
Mesh Leb. const. constant                1.082392200292394                       
Card. mesh Leb. const.                   16                                      
number of interpolation points           8                                       
Mesh interp. nodes degree                2                                       
Mesh interp. nodes factor                2                                       
Card mesh for extracting the nodes       8                                       
================================================================================

Remark: All figures are save in the folder figures

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


## Authors:
    - D. J. Kenne (Code Implementation) <dimitri.kenne@doctoral.uj.edu.pl>
    - Alvise Sommariva <alvise@math.unipd.it>
    - Marco Vianello  <marcov@math.unipd.it>

## Thank you to 
Prof. Dr. hab. Leokadia Biales-Ciez for her help in developing this work.

## References
1. L. Białas-Cież, D. J. Kenne, A. Sommariva, and M. Vianello, ``Chebyshev admissible meshes and Lebesgue constants of complex polynomial projections," J. Comput. Appl. Math., vol. 443, p. 115766, 2024.

