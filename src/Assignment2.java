
public class Assignment2 {

	/*-----------------------
	 *| Part A - tasks 1-11 |
	 * ----------------------*/

	// task 1
	public static boolean isSquareMatrix(boolean[][] matrix) {
		if (matrix == null) {
			return false;
		}
		if (matrix.length == 0) {
			return false;
		}
		for (int row = 0; row < matrix.length; row = row + 1)
			if (matrix.length != matrix[row].length) {
				return false;
			}
		return true;
	}

	// task 2
	public static boolean isSymmetricMatrix(boolean[][] matrix) {
		for (int i = 0; i < matrix.length; i = i + 1)
			for (int j = 0; j < matrix.length; j = j + 1) {
				if (matrix[i][j] != matrix[j][i]) {
					return false;
				}
			}
		return true;
	}

	// task 3
	public static boolean isAntiReflexiveMatrix(boolean[][] matrix) {
		for (int i = 0; i < matrix.length; i = i + 1) {
			if (matrix[i][i] != false)
				return false;
		}
		return true;
	}

	// task 4
	public static boolean isLegalInstance(boolean[][] matrix) {
		// according to the first definition the matrix has to be all of the mentioned below , but first it has to be square.
		if (isSquareMatrix(matrix) == true
				&& (isSymmetricMatrix(matrix) == true & (isAntiReflexiveMatrix(matrix) == true))) {
			return true;
		}
		return false;
	}

	// task 5
	public static boolean isPermutation(int[] array) {
		for (int value = 0; value < array.length; value = value + 1) {
			// checking if we have a value that is bigger than allowed
			if (array[value] > array.length - 1)
				return false;
		}
		for (int i = 0; i < array.length - 1; i = i + 1) {
			for (int j = i + 1; j < array.length; j = j + 1)
				// after we have made sure the values are not too big , we make sure we don't have a certain number twice.
				// if there are no values bigger than allowed and no doubles we exactly have a permutation.
				if (array[i] == array[j])
					return false;
		}
		return true;
	}

	// task 6
	public static boolean hasLegalSteps(boolean[][] flights, int[] tour) {
		// first we make sure the last city has a flight back to city 0.
		if (flights[tour[0]][tour[tour.length - 1]] != true) {
			return false;
		}
		for (int i = 0; i < tour.length - 1; i = i + 1) {
			// making sure that according to tour that stands for out movement between cities , we indeed
			// have a flight between those cities.
			if (flights[tour[i]][tour[i + 1]] != true)
				return false;
		}
		return true;
	}

	// task 7
	//according to the second definition , a solution has to be all of the mentioned below.
	public static boolean isSolution(boolean[][] flights, int[] tour) {
		if (tour.length != flights[0].length)
			throw new IllegalArgumentException("input tour doesn't have a valid length");
		if (isPermutation(tour) == false) {
			return false;
		}
		if (tour[0] != 0) {
			return false;
		}
		if (hasLegalSteps(flights, tour) == false) {
			return false;
		}
		return true;
	}

	// task 8
	public static boolean evaluateLiteral(int literal, boolean[] assign) {
		// first we evaluate the literals.
		boolean LiteralValue;
		if (literal > 0)
			LiteralValue = assign[literal];
		else
			LiteralValue = !assign[-literal];
		return LiteralValue;
	}
	public static boolean evaluateClause(int[] Clause, boolean[] assign) {
		// second , we evaluate each clause
		boolean ClauseValue = false;
		for (int LiteralIndex = 0; LiteralIndex < Clause.length & !ClauseValue; LiteralIndex = LiteralIndex + 1) {
			int literal = Clause[LiteralIndex];
			boolean LiteralValue = evaluateLiteral(literal, assign);
			ClauseValue = LiteralValue;
		}
		return ClauseValue;
	}
	public static boolean evaluate(int[][] cnf, boolean[] assign) {
		// third , we evaluate the whole cnf according to our previous evaluations to see if the assignment fits the cnf.
		boolean cnfValue = true;
		for (int ClauseIndex = 0; ClauseIndex < cnf.length & cnfValue; ClauseIndex = ClauseIndex + 1) {
			cnfValue = evaluateClause(cnf[ClauseIndex], assign);
		}
		return cnfValue;
	}

	

	

	// task 9
	public static int[][] atLeastOne(int[] lits) {
		int[][] cnf = new int[1][lits.length];
		// if we insert all the literals into one clause together in the cnf , at least one of them has to be true.
		cnf[0] = lits;
		return cnf;
	}

	// task 10
	public static int[][] atMostOne(int[] lits) {
		int numOfLits = lits.length;
		int numOfPairs = numOfLits * (numOfLits - 1) / 2;
		int currIndex = 0;
		int[][] cnf = new int[numOfPairs][2];
		for (int i = 0; i < lits.length - 1; i = i + 1) {
			for (int j = i + 1; j < lits.length; j = j + 1, currIndex = currIndex + 1) {
				int[] clause = { -lits[i], -lits[j] };
				cnf[currIndex] = clause;
			}
		}
		return cnf;
	}

	// task 11
	public static int[][] exactlyOne(int[] lits) {
		// combining both atLeastOne one and atMostOne together , and creating a cnf consisting of both the clauses in 
		// them , giving us exactly one.
		int numOfLits = lits.length;
		int numOfPairs = numOfLits * (numOfLits - 1) / 2;
		int numOfClauses = numOfPairs + 1;
		int[][] cnf = new int[numOfClauses][];
		for (int i = 0; i < cnf.length - 1; i = i + 1) {
			cnf[i] = atMostOne(lits)[i];
		}
		cnf[cnf.length - 1] = atLeastOne(lits)[0];
		return cnf;
	}

	/*------------------------
	 *| Part B - tasks 12-20 |
	 * -----------------------*/

	// task 12a
	public static int map(int i, int j, int n) {
		int ans = 0;
		ans = (n * i) + j + 1;
		return ans;
	}

	// task 12b
	public static int[] reverseMap(int k, int n) {
		int i = (k - 1) / n;
		int j = (k - 1) % n;
		int[] ans = { i, j };
		return ans;
	}

	// task 13
	public static int[][] oneCityInEachStep(int n) {
		int[][] matrix = new int [n][n];
		for (int i=0 ; i<n ; i=i+1) {
			for (int j=0;j<n;j=j+1) {
				//creating the map
				matrix[i][j] = map(i,j,n);
			}
		}
		// cnf length has to be exactlyone multiplied by number of steps (rows)
		int[][] cnf = new int [exactlyOne(matrix[0]).length*n][];  
		for (int p = 0 ; p < matrix.length ; p = p+1) {
			// for each row in the map (represents a step) , we make sure exactly one city is visited.
			// then we add all clauses into one cnf , which we declared with the right size.
			int[][] cnfpart = exactlyOne(matrix[p]);
			int startIndex = p*cnfpart.length;
			for (int k=0 ; k<cnfpart.length; k=k+1 , startIndex = startIndex+1) {
				cnf[startIndex] = cnfpart[k];
			}
		}
		return cnf;
	}

	// task 14
	//repeating the proccess in task 13 , only now we make sure each city is visited once.
	public static int[][] eachCityIsVisitedOnce(int n) {
		int[][] matrix = new int [n][n];
		for (int j=0 ; j<n ; j=j+1) {
			for (int i=0;i<n;i=i+1) {
				matrix[j][i] = map(i,j,n);
			}
		}
		int[][] cnf = new int [exactlyOne(matrix[0]).length*n][];
		for (int p = 0 ; p < matrix.length ; p = p+1) {
			int[][] cnfpart = exactlyOne(matrix[p]);
			int startIndex = p*cnfpart.length;
			for (int k=0 ; k<cnfpart.length; k=k+1 , startIndex = startIndex+1) {
				cnf[startIndex] = cnfpart[k];
			}
		}
		return cnf;
	}

	// task 15
	public static int[][] fixSourceCity(int n) {
		int[][] cnf = {{map(0,0,n)}};
		return cnf;
	}

	// task 16
	public static int[][] noIllegalSteps(boolean[][] flights) {
		int n=flights.length;
		int count = 0;
		for (int i=0;i<n-1;i=i+1) {
			for (int j=i+1; j<n;j=j+1) {
				//checking how many false values we have above all ID values - (1,1) (2,2) etc , because
				// the matrix has to be symmetric , and false values in ID values don't count because they represent
				//a flight from a city to itself.
				if (flights[i][j] == false)
					count=count+2;
			}
		}
		int [][] finalcnf = new int [count*n][2];
		int cnfindex = 0;
		for (int i = 0 ; i<n ; i=i+1) {
			for (int j = 0; j<n ; j=j+1) {
				for (int k = 0 ; k<n ; k=k+1) {
					// If we have found a false value that isn't an ID value , we make sure we don't make a step
					// between those cities. after creating the clauses that prevent those steps , we concatenate them into
					// one cnf formula
					if (flights[j][k] == false & (j!=k)) {
						int[] clause = {-(map(i,j,n)) , -(map((i+1)%n,k,n))};
						finalcnf[cnfindex] = clause;
						cnfindex = cnfindex + 1;
					}
				}
			}
		}
		return finalcnf;
	}	
	// task 17
	// I created a function that concatenates two cnf formulas to ease on the process and make the code readable.
	public static int [][] CNFcombiner (int[][] cnf1 , int[][] cnf2) {
		int [][] combinedCNF = new int [cnf1.length + cnf2.length][];
		for (int i = 0 ; i<cnf1.length ; i=i+1) {
			combinedCNF[i] = cnf1[i];
		}
		int currIndex = cnf1.length;
		for (int i = 0 ; i<cnf2.length ; i=i+1 , currIndex=currIndex+1) {
			combinedCNF[currIndex] = cnf2[i];
		}
		return combinedCNF;
	}
	public static int[][] encode(boolean[][] flights) {
		int n = flights.length;
		int [][] cnf1 = oneCityInEachStep (n);
		int [][] cnf2 = eachCityIsVisitedOnce (n);
		int [][] cnf3 = fixSourceCity (n);
		int [][] cnf4 = noIllegalSteps (flights);
		// using the cnf combining function to concatenate all previous terms
		int [][] finalcnf = CNFcombiner(cnf1,cnf2);
		finalcnf = CNFcombiner(finalcnf,cnf3);
		finalcnf = CNFcombiner(finalcnf,cnf4);
		return finalcnf;
	}

	// task 18
	public static int[] decode(boolean[] assignment, int n) {
		if (assignment == null | assignment.length != n*n+1)
			throw new IllegalArgumentException();
		int[] tour = new int [n];
		int tourIndex=0;
		for (int k=0; k<assignment.length ; k=k+1) {
			// if we found a k that is true , we extract j that stands for the city.
			if (assignment[k] == true) {
			int j = (k-1)%n;
			tour[tourIndex] = j;
			tourIndex = tourIndex + 1;
			}
		}
		return tour;
	}

	// task19
	public static int[] solve(boolean[][] flights) {
		if (isLegalInstance(flights) == false)
			throw new IllegalArgumentException();
		int n = flights.length;
		int nVars = n*n;
		SATSolver.init(nVars);
		int[][] cnf = encode(flights);
		SATSolver.addClauses(cnf);
		boolean[] assignment = SATSolver.getSolution();
		if (assignment == null)
			throw new IllegalArgumentException();
		if (assignment.length == nVars+1) {
			int[] tour = decode(assignment,n);
			if (isSolution(flights,tour) == true) {
				return tour;
			}
			else throw new IllegalArgumentException();
		}
		else return null;	
	}

	// task20
	public static boolean solve2(boolean[][] flights) {
		boolean ans = true;
		if (isLegalInstance(flights) == false)
			throw new IllegalArgumentException();
		int[] tour = solve(flights);
		if (tour == null) 
			ans = false;
		int n = flights.length;
		//after we had found a solution , we want to make sure the equal solution is not used so I created the equal tour.
		int[] reversedTour = new int [tour.length];
		for (int i = 1 ; i<tour.length ; i=i+1) {
			reversedTour[i] = tour[tour.length-i];
		}
		//after we have both tours , we want to negate them so they are not used the next time we try to find
		//a different path to use.
		int[] clause1 = new int [tour.length-1];
		int[] clause2 = new int [reversedTour.length-1];
		for (int i = 1 ; i<tour.length ; i = i+1) {
			clause1[i-1] = -map(i,tour[i],n);
			clause2[i-1] = -map(i,reversedTour[i],n);
		}
		int nVars = n*n;
		SATSolver.init(nVars);
		//after we have both the negations of the first solution and the equal tour
		//we add both clauses , and then add our first cnf.
		SATSolver.addClause(clause1);
		SATSolver.addClause(clause2);
		int[][] cnf = encode(flights);
		SATSolver.addClauses(cnf);
		boolean[] assignment = SATSolver.getSolution();
		if (assignment == null)
			throw new IllegalArgumentException();
		if (assignment.length == nVars+1) {
			int[] secondtour = decode(assignment,n);
			if (isSolution(flights,secondtour) == true) {
				ans = true;
			}
			else ans = false;
		}
		else ans = false;
		return ans;
	}
}