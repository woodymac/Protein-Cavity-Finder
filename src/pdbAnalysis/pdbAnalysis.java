package pdbAnalysis;

// This block of code provides the necessary imports for the pdbAnalysis file.
import java.util.ArrayList;
import java.util.Scanner;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.lang.Double;

/**
 * Written in Eclipse.
 * 
 * John Heberle, Woody McIlhenny, Zach King, Matthew Schaefer
 * 
 * This class reads the input file and processes the data. The data processing
 * is
 * handled by storing each atom read in by the scanner into an array of atoms.
 * This array is then turned into a three dimensional array of coordinates which
 * is then passed through to other processing classes.
 */
public class pdbAnalysis {

    /**
     * 
     * @param args
     *            [input file name][output file name][radius][resolution]
     * 
     * @throws FileNotFoundException
     *             If the input is not in the directory
     * @throws UnsupportedEncodingException 
     * 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {

        Atom[] atomIn = new Atom[5000];
        String pdbFileName;
        String outputFile;
        double probeR;
        double rez;
        ArrayList<Atom> atomArray = new ArrayList<Atom>();

        // If there are no defined arguments, ask user directly for input.
        if (args.length == 0) {
            System.out.println(
                "Input the input file, output file, probe size, and resolution:");
            System.out.println("(Ex: input.pdb output.pdb 0.25 0.5)");

            Scanner pdbInScan = new Scanner(System.in);
            pdbFileName = pdbInScan.next();
            outputFile = pdbInScan.next();
            String radius = pdbInScan.next();
            probeR = Double.parseDouble(radius);
            String resolution = pdbInScan.next();
            rez = Double.parseDouble(resolution);
        }
        else {
            pdbFileName = args[0];
            outputFile = args[1];
            probeR = Double.parseDouble(args[2]);
            rez = Double.parseDouble(args[3]);
        }

        // The following code sets the primary data acquisition variables and
        // the scanning functionality.

        File pdbFile = new File(pdbFileName);

        // The following code sets the initial variables used to read in the
        // information from the input file.
        Scanner pdbAtoms = new Scanner(pdbFile);
        String atomInfoLine;
        String[] atomInfo;
        int atomNum = 0;
        double minValX = 0;
        double minValY = 0;
        double minValZ = 0;
        double maxValX = 0;
        double maxValY = 0;
        double maxValZ = 0;

        // The follwoing block of code reads through the input file and and
        // stores
        // the atoms as an atom in the atomIn array.
        while (pdbAtoms.hasNext()) {

            // This initializes the variables used to store the location of the
            // atom.
            double x = 0;
            double y = 0;
            double z = 0;
            double r = 0;

            // This reads in the next line from the pdb file and extracts the
            // relevant information.
            atomInfoLine = pdbAtoms.nextLine();
            String atomInfoLineP = atomInfoLine.substring(33, atomInfoLine
                .length());
            Scanner line = new Scanner(atomInfoLineP);

            x = line.nextDouble();
            if (x < minValX) {
                minValX = x;
            }
            if (x > maxValX) {
                maxValX = x;
            }
            y = line.nextDouble();
            if (y < minValY) {
                minValY = y;
            }
            if (y > maxValY) {
                maxValY = y;
            }
            z = line.nextDouble();
            if (z < minValZ) {
                minValZ = z;
            }
            if (z > maxValZ) {
                maxValZ = z;
            }
            line.nextDouble();
            r = line.nextDouble();

            Atom a = new Atom(x, y, z, r);
            atomArray.add(a);
            atomIn[atomNum] = a;
            atomNum++;
        }

        /**
         * Finds x, y, and z ranges of the gridSpace
         **/
        int xSize = ((int)((maxValX - minValX)/rez) + 1);
        int ySize = ((int)((maxValY - minValY)/rez) + 1);
        int zSize = ((int)((maxValZ - minValZ)/rez) + 1);

        /**
         * Creates new gridSpace
         **/
        int[][][] gridSpace = new int[zSize][ySize][xSize];

        /**
         * Loops through each atom, finding its z-range
         **/
        for (int i = 0; i < atomNum; i++) {
            Atom searchAtom = atomIn[i];
            double lowestX = minValX;
            double lowestY = minValY;
            double lowestZ = minValZ;
            double x = searchAtom.x;
            double y = searchAtom.y;
            double z = searchAtom.z;
            double r = searchAtom.r;

            int topZLayer = (int)((z - r - lowestZ) / rez);
            int midZLayer = (int)((z - lowestZ) / rez);
            int bottomZLayer = (int)((z + r - lowestZ) / rez);
            if (topZLayer < 0) {
                topZLayer = 0;
            }

            int currentLayer = topZLayer;

            /**
             * Loops through each 2D slice that the atom touches
             */
            while (currentLayer < bottomZLayer && currentLayer < (int)((maxValZ
                - minValZ) * 1/rez)) {

                /**
                 * Determines the radius of the atom at the particular slice.
                 */
                double distance = Math.abs(z - (currentLayer * rez - Math.abs(
                    lowestZ)));
                double rNew = r * Math.sqrt(1 - (distance / r) * (distance
                    / r));

                /**
                 * The lowest x-value that the atom touches at the slice, or 0
                 * if the atom exceeds the gridSpace.
                 */
                int xLow = (int)(((x - rNew - lowestX) / rez));
                if (xLow < 0) {
                    xLow = 0;
                }
                int xHigh = (int)(((x + rNew - lowestX) / rez));
                if (xHigh > (int)(maxValX/rez)) {
                    xHigh = (int)(maxValX/rez);
                }

                /**
                 * The lowest y-value that the atom touches at the slice, or 0
                 * if the atom exceeds the gridSpace.
                 */
                int yLow = (int)(((y - rNew - lowestY) / rez));
                if (yLow < 0) {
                    yLow = 0;
                }
                int yHigh = (int)(((y + rNew - lowestY) / rez));
                if (yHigh > (int)(maxValY/rez)) {
                    yHigh = (int)(maxValY/rez);
                }

                /**
                 * Determines the voxels that are touched by the atom on a
                 * particular 2D slice.
                 */
                for (int gridY = yLow; gridY < yHigh; gridY++) {
                    for (int gridX = xLow; gridX < xHigh; gridX++) {
                        gridSpace[currentLayer][gridY][gridX] = 1;
                    }
                }

                currentLayer++;
            }
        }

        /**
         * Finds the cavities in the empty spaces.
         */
        CavityFinder finder = new CavityFinder(gridSpace);
        int[][] cavities = finder.getCavities();

        /**
         * Adds the cavities to the existing gridSpace, where they are
         * represented as -1's.
         */
        for (int c = 0; c < cavities.length; c++) {
            int gridX = cavities[c][0];
            int gridY = cavities[c][1];
            int gridZ = cavities[c][2];
            gridSpace[gridZ][gridY][gridX] = -1;
        }

        double lowestX = minValX;
        double lowestY = minValY;
        double lowestZ = minValZ;
        int probeCount = 0;
        boolean empty = true;
        int counter = 0;
        ArrayList<Integer> finalArray = new ArrayList<Integer>();

        /**
         * Iterates through the list of cavities, testing to see if a probe
         * sphere of the specified radius can be centered in one.
         */
        for (int i = 0; i < cavities.length; i++) {
            counter = 0;
            int x = cavities[i][0];
            int y = cavities[i][1];
            int z = cavities[i][2];

            /**
             * These int's represent the highest and lowest x, y, and z
             * positions
             * in the 3D-array that the probe sphere would touch. If the probe
             * sphere would extend outside of the gridspace, the maximum or
             * minimum
             * value in that direction is set to the final gridspace position.
             */
            int topZLayer = (int)((z * rez - probeR - lowestZ) / rez);
            if (topZLayer < 0) {
                topZLayer = 0;
            }
            int bottomZLayer = (int)((z * rez + probeR - lowestZ) / rez);
            if (bottomZLayer > maxValZ/rez) {
                bottomZLayer = (int)(maxValZ/rez);
            }
            int zRange = bottomZLayer - topZLayer;
            int topXLayer = (int)((x * rez - probeR - lowestX) / rez);
            if (topXLayer < 0) {
                topXLayer = 0;
            }
            int bottomXLayer = (int)((x * rez + probeR - lowestX) / rez);
            if (bottomXLayer > maxValX/rez) {
                bottomXLayer = (int)(maxValX/rez);
            }
            int xRange = bottomXLayer - topXLayer;
            int topYLayer = (int)((y * rez - probeR - lowestY) / rez);
            if (topYLayer < 0) {
                topYLayer = 0;
            }
            int bottomYLayer = (int)((y * rez + probeR - lowestY) / rez);
            if (bottomYLayer > maxValY/rez) {
                bottomYLayer = (int)(maxValY/rez);
            }
            int yRange = bottomYLayer - topYLayer;
            int[][] resultArray = new int[xRange * yRange * zRange][3];
            empty = true;

            /**
             * Searches the area around the potential probe sphere for
             * filled voxels. If one is found, the probe sphere is discarded
             * and the next cavity is tried.
             */
            for (int zLayer = topZLayer; zLayer < bottomZLayer - 1; zLayer++) {
                for (int yLayer = topYLayer; yLayer < bottomYLayer
                    - 1; yLayer++) {
                    for (int xLayer = topXLayer; xLayer < bottomXLayer
                        - 1; xLayer++) {
                        if (gridSpace[zLayer][yLayer][xLayer] == -1) {
                            resultArray[counter][0] = xLayer;
                            resultArray[counter][1] = yLayer;
                            resultArray[counter][2] = zLayer;
                            counter++;
                        }
                        else {
                            empty = false;
                        }
                    }
                }
            }

            /**
             * If the area around a potential probe sphere was found to be
             * empty, the number of probe spheres is incremented, the gridspace
             * is updated to reflect the now-probe-sphere-filled positions, and
             * the probe sphere's center coordinates are added to an arrayList.
             */
            if (empty) {
                probeCount++;
                for (int res = 0; res < resultArray.length; res++) {
                    gridSpace[resultArray[res][2]][resultArray[res][1]][resultArray[res][0]] =
                        -2;
                }
                finalArray.add(x);
                finalArray.add(y);
                finalArray.add(z);
            }

        }

        /**
         * The arrayList of probe spheres, formed from the coordinates
         * found above.
         */
        ArrayList<Atom> probeArray = new ArrayList<Atom>();
        int atomCount = 0;
        while (atomCount < finalArray.size()) {
            double newX = finalArray.get(atomCount) * rez;
            double newY = finalArray.get(atomCount + 1) * rez;
            double newZ = finalArray.get(atomCount + 2) * rez;
            Atom newAtom = new Atom(newX, newY, newZ, probeR);
            probeArray.add(newAtom);
            atomCount = atomCount + 3;

        }

        /**
         * Creates a pdb file of a specified name that contains
         * the list of existing atoms and the probe spheres.
         */
        File outputPDBFile = new File(outputFile);
        PrintWriter writer = new PrintWriter(outputPDBFile, "UTF-8");
        for (int i = 0; i < atomArray.size(); i++) {
            StringBuilder number = new StringBuilder(Integer.toString(i + 1));
            while (number.length() < 7) {
                number.insert(0, " ");
            }
            StringBuilder atomX = new StringBuilder(Double.toString(atomArray
                .get(i).x));
            atomX.append("0");
            atomX.append("0");
            StringBuilder atomY = new StringBuilder(Double.toString(atomArray
                .get(i).y));
            atomY.append("0");
            atomY.append("0");
            StringBuilder atomZ = new StringBuilder(Double.toString(atomArray
                .get(i).z));
            atomZ.append("0");
            atomZ.append("0");
            StringBuilder atomR = new StringBuilder(Double.toString(atomArray
                .get(i).r));
            atomR.append("0");
            writer.println("ATOM" + number + " SPH  ION" + number + "     "
                + atomX + "  " + atomY + "  " + atomZ + " 0.00  " + atomR);
        }
        for (int i = 0; i < probeArray.size(); i++) {
            StringBuilder number = new StringBuilder(Integer.toString(i
                + 10001));
            while (number.length() < 7) {
                number.insert(0, " ");
            }
            StringBuilder atomX = new StringBuilder(Double.toString(probeArray
                .get(i).x));
            atomX.append("0");
            atomX.append("0");

            StringBuilder atomY = new StringBuilder(Double.toString(probeArray
                .get(i).y));
            atomY.append("0");
            atomY.append("0");

            StringBuilder atomZ = new StringBuilder(Double.toString(probeArray
                .get(i).z));
            atomZ.append("0");
            atomZ.append("0");

            StringBuilder atomR = new StringBuilder(Double.toString(probeArray
                .get(i).r));
            atomR.append("0");
            atomR.append("0");

            writer.println("ATOM" + number + "  MC  CAV" + number + "     "
                + atomX + "  " + atomY + "  " + atomZ + " 0.00  " + atomR);

        }

        writer.close();
        System.out.println();
        System.out.println("Number of Probes: " + probeCount);

    }


    /**
     * The Atom class contains an atom's coordinates and radius.
     */
    public static class Atom {
        public double x;
        public double y;
        public double z;
        public double r;


        Atom(double x, double y, double z, double r) {
            this.x = x;
            this.y = y;
            this.z = z;
            this.r = r;
        }


        /**
         * Returns the x coordinate.
         */
        public double getX() {
            return x;
        }


        /**
         * Returns the y coordinate.
         */
        public double getY() {
            return y;
        }


        /**
         * Returns the z coordinate.
         */
        public double getZ() {
            return z;
        }


        /**
         * Returns the atom's coordinates.
         */
        public String toString() {
            return x + " " + y + " " + z;
        }
    }


    /**
     * Takes a 3d matrix representing the protein, and returns an array
     * containing
     * the
     * coordinates of the cavities.
     * 
     * @author Woody McIlhenny
     *
     */
    public static class CavityFinder {

        // Matrix representing protein, [z][y][x]
        int[][][] grid;

        // List of void spaces at extreme x, y, and z
        ArrayList<Atom> boundList;

        // List of void spaces identified as NOT being cavities
        ArrayList<Atom> voidList;

        // Once getCavitites is run, will only contain cavities
        ArrayList<Atom> cavityList;


        /**
         * Creates a new CavityFinder object.
         * 
         * @param 3D
         *            matrix representing the protein. Each array element
         *            represents
         *            a coordinate, formatted as [z][y][x]. A 1 is an atom at
         *            those
         *            coordinates, a 0 is a null space.
         */
        public CavityFinder(int[][][] gridSpace) {
            grid = gridSpace;
            boundList = new ArrayList<Atom>();
            voidList = new ArrayList<Atom>();
            cavityList = new ArrayList<Atom>();

            // Print out a visual representation of the protein
            System.out.println("Protein//////////////////////");
            for (int z = 0; z < grid.length; z++) {
                for (int y = 0; y < grid[z].length; y++) {
                    for (int x = 0; x < grid[z][y].length; x++) {
                        System.out.print(grid[z][y][x]);
                    }
                    System.out.println();
                }
                System.out.println();
            }
            System.out.println("//////////////////////////////");
        }


        /**
         * Returns an array of the coordinates for all cavities in the following
         * format:
         * [x1][y1][z1]
         * [x2][y2][z2]
         * [x3][y3][z3]...
         */
        public int[][] getCavities() {
            populateBoundList();
            populateCavityList();
            checkIfCavityFront();
            checkIfCavityBack();

            System.out.println("Cavities:");
            // make array
            int[][] cavityArray = new int[cavityList.size()][3];
            for (int i = 0; i < cavityList.size(); i++) {
                cavityArray[i][0] = (int)cavityList.get(i).getX();
                cavityArray[i][1] = (int)cavityList.get(i).getY();
                cavityArray[i][2] = (int)cavityList.get(i).getZ();
                System.out.print(i + 1 + ". x=" + cavityArray[i][0] + " ");
                System.out.print("y=" + cavityArray[i][1] + " ");
                System.out.print("z=" + cavityArray[i][2] + " ");
                System.out.println();
            }
            if (cavityArray.length == 0) {
                System.out.println("No Cavitities");
            }
            return cavityArray;
        }


        /**
         * Creates list of voids on x, y, and z extremes
         */
        private void populateBoundList() {

            // add z-extremes
            int z = 0;
            for (int y = 0; y < grid[z].length; y++) {
                for (int x = 0; x < grid[z][y].length; x++) {
                    if (grid[z][y][x] == 0 && !preexistingAtom(boundList, x, y,
                        z)) {
                        boundList.add(new Atom(x, y, z, 0));
                    }
                }
            }
            z = grid.length - 1;
            for (int y = 0; y < grid[z].length; y++) {
                for (int x = 0; x < grid[z][y].length; x++) {
                    if (grid[z][y][x] == 0 && !preexistingAtom(boundList, x, y,
                        z)) {
                        boundList.add(new Atom(x, y, z, 0));
                    }
                }
            }

            // add y-extremes
            int y = 0;
            for (z = 0; z < grid.length; z++) {
                for (int x = 0; x < grid[z][y].length; x++) {
                    if (grid[z][y][x] == 0 && !preexistingAtom(boundList, x, y,
                        z)) {
                        boundList.add(new Atom(x, y, z, 0));
                    }
                }
            }
            for (z = 0; z < grid.length; z++) {
                y = grid[z].length - 1;
                for (int x = 0; x < grid[z][y].length; x++) {
                    if (grid[z][y][x] == 0 && !preexistingAtom(boundList, x, y,
                        z)) {
                        boundList.add(new Atom(x, y, z, 0));
                    }
                }
            }

            // add x-extremes
            int x = 0;
            for (z = 0; z < grid.length; z++) {
                for (y = 0; y < grid[z].length; y++) {
                    if (grid[z][y][x] == 0) {
                        boundList.add(new Atom(x, y, z, 0));
                    }
                }
            }
            for (z = 0; z < grid.length; z++) {
                for (y = 0; y < grid[z].length; y++) {
                    x = grid[z][y].length - 1;
                    if (grid[z][y][x] == 0 && !preexistingAtom(boundList, x, y,
                        z)) {
                        boundList.add(new Atom(x, y, z, 0));
                    }
                }
            }
        }


        /**
         * Determines if any of the atoms in the given list contains
         * an atom with the given coordinates.
         */
        private boolean preexistingAtom(
            ArrayList<Atom> list,
            int x,
            int y,
            int z) {
            for (Atom atom : list) {
                if (atom.getX() == x && atom.getY() == y && atom.getZ() == z) {
                    return true;
                }
            }
            return false;
        }


        /**
         * Creates the list of void spaces, which are zeroes in the
         * grid array. They are considered cavities initially, until found to be
         * otherwise.
         */
        private void populateCavityList() {
            // assumes grid is at least 3x3x3
            for (int z = 1; z < grid.length - 1; z++) {
                for (int y = 1; y < grid[z].length - 1; y++) {
                    for (int x = 1; x < grid[z][y].length - 1; x++) {
                        if (grid[z][y][x] == 0) {
                            cavityList.add(new Atom(x, y, z, 0));
                        }
                    }
                }
            }
        }


        /**
         * Checks through the voidList to see if void space is connected to
         * outside.
         * Starts at front of list.
         */
        private void checkIfCavityFront() {
            for (int i = 0; i < cavityList.size(); i++) {
                Atom atom = cavityList.get(i);
                int thisX = (int)atom.getX();
                int thisY = (int)atom.getY();
                int thisZ = (int)atom.getZ();

                // check if surrounding voids are on the boundary
                if (preexistingAtom(boundList, thisX + 1, thisY, thisZ)
                    || preexistingAtom(boundList, thisX, thisY + 1, thisZ)
                    || preexistingAtom(boundList, thisX, thisY, thisZ + 1)
                    || preexistingAtom(boundList, thisX - 1, thisY, thisZ)
                    || preexistingAtom(boundList, thisX, thisY - 1, thisZ)
                    || preexistingAtom(boundList, thisX, thisY, thisZ - 1)) {
                    voidList.add(atom);
                    cavityList.remove(atom);
                    i--;
                }
                // check if surrounding voids have been deemed not cavities
                else if (preexistingAtom(voidList, thisX + 1, thisY, thisZ)
                    || preexistingAtom(voidList, thisX, thisY + 1, thisZ)
                    || preexistingAtom(voidList, thisX, thisY, thisZ + 1)
                    || preexistingAtom(voidList, thisX - 1, thisY, thisZ)
                    || preexistingAtom(voidList, thisX, thisY - 1, thisZ)
                    || preexistingAtom(voidList, thisX, thisY, thisZ - 1)) {
                    voidList.add(atom);
                    cavityList.remove(atom);
                    i--;
                }
            }
        }


        /**
         * Checks through the voidList to see if void space is connected to
         * outside.
         * Starts at back of list.
         */
        private void checkIfCavityBack() {
            for (int i = cavityList.size() - 1; i >= 0; i--) {
                Atom atom = cavityList.get(i);
                int thisX = (int)atom.getX();
                int thisY = (int)atom.getY();
                int thisZ = (int)atom.getZ();

                // check if surrounding voids are on the boundary
                if (preexistingAtom(boundList, thisX + 1, thisY, thisZ)
                    || preexistingAtom(boundList, thisX, thisY + 1, thisZ)
                    || preexistingAtom(boundList, thisX, thisY, thisZ + 1)
                    || preexistingAtom(boundList, thisX - 1, thisY, thisZ)
                    || preexistingAtom(boundList, thisX, thisY - 1, thisZ)
                    || preexistingAtom(boundList, thisX, thisY, thisZ - 1)) {
                    voidList.add(atom);
                    cavityList.remove(atom);
                }
                // check if surrounding voids have been deemed not cavities
                else if (preexistingAtom(voidList, thisX + 1, thisY, thisZ)
                    || preexistingAtom(voidList, thisX, thisY + 1, thisZ)
                    || preexistingAtom(voidList, thisX, thisY, thisZ + 1)
                    || preexistingAtom(voidList, thisX - 1, thisY, thisZ)
                    || preexistingAtom(voidList, thisX, thisY - 1, thisZ)
                    || preexistingAtom(voidList, thisX, thisY, thisZ - 1)) {
                    voidList.add(atom);
                    cavityList.remove(atom);
                }
            }
        }
    }

}
