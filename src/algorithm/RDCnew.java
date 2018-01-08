//
package algorithm;

import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import coser.common.SimpleTool;

import tree.Tree;
import weka.clusterers.EM;
import weka.clusterers.FarthestFirst;
import weka.clusterers.HierarchicalClusterer;
import weka.clusterers.SimpleKMeans;
import weka.core.Instance;
import weka.core.Instances;
import weka.filters.unsupervised.attribute.Remove;

public class RDCnew extends Instances {
	private static final long serialVersionUID = -82397417233430226L;
	private static String DataSet = "sonar";
	private static String fileAddress = "/Users/dengsiyu/eclipse/workspace/Coser_Triple/data/wdbc.arff";
	public static int numClass;
	///Users/dengsiyu/eclipse/workspace/Coser_Triple/data

	protected double[][] similarityMatrix;

	public boolean[] alreadyClassified;

	protected double[][] distanceMatrix;
	protected double[][] powerMatrix;
	protected double[][] modifiedMatrix;
	protected double[] rankVector;
	protected double zuni = 0.95;
	protected double threshold = 0.1;
	protected static double percentage = 0.4;
	protected int times = 0;
	double treethreshold = 0.001;
	protected double[] thresholdVector;
	protected int[] neighborCountVector;
	protected int[][] neighborMatrix;

	public int numVote;
	public int numPredict;
	public static int numTeach;

	public int[] predictLabel;

	/**
	 * The block information.
	 */
	int[][] blockInformation;
	public static int[] clusterNumbers;
	public int clusterNum = 0;

	public RDCnew(Reader paraReader) throws Exception, IOException {
		super(paraReader);
	}// Of constructor

	public RDCnew(Instances instances) throws Exception, IOException {
		super(instances);
	}// Of constructor

	public static double getSimilarityValue(Instance former, Instance latter) {
		int count = 0;
		for (int i = 0; i < former.numAttributes() - 1; i++) {
			int valueFormer = (int) former.value(i);
			int valueLatter = (int) latter.value(i);
			if (valueFormer == valueLatter) {
				count++;
			} // Of if
		} // Of for i

		return (count + 0.0) / (former.numAttributes() - 1);
	}// Of getSimilarityValue

	public static double getDistanceValue(Instance former, Instance latter) {

		double tempDistance = 0;

		for (int i = 0; i < former.numAttributes() - 1; i++) {
			double a = former.value(i);
			double b = latter.value(i);
			tempDistance += Math.abs(a - b);
		} // Of for i

		return tempDistance;
	}// Of getDistcaneValue

	public int getMaxIndex(int[] paraArray) {
		int maxIndex = 0;
		int tempIndex = 0;
		int max = paraArray[0];

		for (int i = 0; i < paraArray.length; i++) {
			if (paraArray[i] > max) {
				max = paraArray[i];
				tempIndex = i;
			} // of if
		} // of for i
		maxIndex = tempIndex;
		return maxIndex;
	}// of getMaxIndex

	public void getSimilarityMatrix() {
		similarityMatrix = new double[numInstances()][numInstances()];
		for (int i = 0; i < numInstances(); i++) {
			for (int j = i; j < numInstances(); j++) {
				similarityMatrix[i][j] = 1.0 / (distanceMatrix[i][j] + 1.0);
				similarityMatrix[j][i] = similarityMatrix[i][j];
			} // Of for j
		} // Of for i
		System.out.println("查看相似度矩阵" + Arrays.deepToString(similarityMatrix));

		// System.out.println("184 simi"+Arrays.toString(similarityMatrix[184]));
	}// Of getSimilarityMatrix

	public void getDistanceMatrix() {
		distanceMatrix = new double[numInstances()][numInstances()];
		for (int i = 0; i < numInstances(); i++) {
			for (int j = i; j < numInstances(); j++) {
				distanceMatrix[i][j] = getDistanceValue(instance(i), instance(j));
				distanceMatrix[j][i] = distanceMatrix[i][j];
			} // Of for j
		} // Of for i

		System.out.println("distance" + Arrays.deepToString(distanceMatrix));
	}// Of getSimilarityMatrix

	public void asignMostSimilar() {
		double[] mostSimilarVector = new double[numInstances()];

		// // 为相似度大于0.8 的实例分配分值
		// boolean[][] statuMatrix = new boolean[numInstances()][numInstances()];
		// int[] countVector = new int[numInstances()];
		// for (int i = 0; i < numInstances(); i++) {
		// for (int j = 0; j < numInstances(); j++) {
		// if (similarityMatrix[i][j] >= 0.88) {
		// countVector[i]++;
		// statuMatrix[i][j] = true;
		// } // Of if
		// } // Of for j
		// } // Of for i
		//
		// powerMatrix = new double[numInstances()][numInstances()];
		// for (int i = 0; i < numInstances(); i++) {
		// double value = (1 + 0.0) / countVector[i];
		// for (int j = 0; j < numInstances(); j++) {
		// if (statuMatrix[i][j]) {
		// powerMatrix[i][j] = value;
		// } // Of if
		// } // Of for j
		// } // Of for i

		for (int i = 0; i < numInstances(); i++) {
			double max = Double.MIN_VALUE;
			double secondMax = Double.MIN_VALUE;
			;
			for (int j = 0; j < numInstances(); j++) {
				if (i != j) {
					if (similarityMatrix[i][j] > secondMax) {
						if (similarityMatrix[i][j] > max) {
							double temp = max;
							max = similarityMatrix[i][j];
							secondMax = temp;
						} else {
							secondMax = similarityMatrix[i][j];
						}
						// max = similarityMatrix[i][j];
					} // Of if
				} // Of if
			} // Of for j
			mostSimilarVector[i] = secondMax;
		} // Of for i

		boolean[][] statuMatrix = new boolean[numInstances()][numInstances()];
		int[] countVector = new int[numInstances()];
		for (int i = 0; i < numInstances(); i++) {
			for (int j = 0; j < numInstances(); j++) {
				if (similarityMatrix[i][j] >= mostSimilarVector[i]) {
					countVector[i]++;
					statuMatrix[i][j] = true;
				} // Of if
			} // Of for j
		} // Of for i

		powerMatrix = new double[numInstances()][numInstances()];
		for (int i = 0; i < numInstances(); i++) {
			double value = (1 + 0.0) / countVector[i];
			for (int j = 0; j < numInstances(); j++) {
				if (statuMatrix[i][j]) {
					powerMatrix[i][j] = value;
				} // Of if
			} // Of for j
		} // Of for i

	}// Of asignMostSimilar

	public void asignMostClose() {
		double[] mostfarVector = new double[numInstances()];
		double[] mostCloserVector = new double[numInstances()];
		// 实例i与其他实例的最远距离
		for (int i = 0; i < numInstances(); i++) {
			double max = Double.MIN_VALUE;
			for (int j = 0; j < numInstances(); j++) {
				if (i != j) {
					if (distanceMatrix[i][j] > max) {
						max = distanceMatrix[i][j];
					} // Of if
				} // Of if
			} // Of for j
			mostfarVector[i] = max;
		} // Of for i
			// 实例i与其他实例的最近距离
		for (int i = 0; i < numInstances(); i++) {
			double min = Double.MAX_VALUE;
			for (int j = 0; j < numInstances(); j++) {
				if (i != j) {
					if (distanceMatrix[i][j] < min) {
						min = distanceMatrix[i][j];
					} // Of if
				} // Of if
			} // Of for j
			mostCloserVector[i] = min;
		} // Of for i
		int[] sortdistanceMatrix = new int[numInstances()];

		boolean[][] statuMatrix = new boolean[numInstances()][numInstances()];
		int[] countVector = new int[numInstances()];
		for (int i = 0; i < numInstances(); i++) {
			sortdistanceMatrix = mergeSortToIndices(distanceMatrix[i]);

			for (int j = 0; j < 5; j++) {

				statuMatrix[i][sortdistanceMatrix[numInstances() - j - 1]] = true;

			} // Of for j
		} // Of for i

		powerMatrix = new double[numInstances()][numInstances()];
		for (int i = 0; i < numInstances(); i++) {
			double value = (1 + 0.0) / 5;
			for (int j = 0; j < numInstances(); j++) {
				if (statuMatrix[i][j]) {
					powerMatrix[i][j] = value;
				} // Of if
			} // Of for j
		} // Of for i

		// System.out.println("184 power"+Arrays.toString(powerMatrix[184]));

		// 为相似度大于0.8 的实例分配分值
		// boolean[][] statuMatrix = new boolean[numInstances()][numInstances()];
		// int[] countVector = new int[numInstances()];
		// for (int i = 0; i < numInstances(); i ++) {
		// for (int j = 0; j < numInstances(); j ++) {
		// if (similarityMatrix[i][j] >= 0.8*max) {
		// countVector[i] ++;
		// statuMatrix[i][j] = true;
		// }// Of if
		// }// Of for j
		// }// Of for i
		//
		// powerMatrix = new double[numInstances()][numInstances()];
		// for (int i = 0; i < numInstances(); i ++) {
		// double value = (1 + 0.0) / countVector[i];
		// for (int j = 0; j < numInstances(); j ++) {
		// if (statuMatrix[i][j]) {
		// powerMatrix[i][j] = value;
		// }// Of if
		// }// Of for j
		// }// Of for i

	}// Of asignMostSimilar

	public int getMax(int[] paracluster) {

		int tempMax = Integer.MIN_VALUE;
		for (int i = 0; i < paracluster.length; i++) {

			if (paracluster[i] > tempMax)

				tempMax = clusterNumbers[i];

		} // of for i
		return tempMax;

	}// of

	public static double getPredictionAccuracy(int[] para1, int[] para2) {
		double tempInCorrect = 0;
		// System.out.println("Incorrectly classified instances:");
		for (int i = 0; i < para1.length; i++) {
			if (para1[i] != para2[i]) {
				tempInCorrect++;
				System.out.print("" + i + ", ");
			} // Of if
		} // Of for i
		System.out.println();
		System.out.println("This is the incorrect:\r\n" + tempInCorrect);
		System.out.println("total" + para1.length);
		
	    System.out.println("精度" + (para1.length - tempInCorrect - numTeach) / (para1.length-numTeach));
		return tempInCorrect / para1.length;

	}// Of getPredictionAccuracy

	public void revertPowerMatrix() {
		double[][] tempMatrix = new double[numInstances()][numInstances()];
		for (int i = 0; i < numInstances(); i++) {
			for (int j = 0; j < numInstances(); j++) {
				tempMatrix[i][j] = powerMatrix[j][i];
			} // Of for j
		} // Of for i
		powerMatrix = tempMatrix;

		System.out.println("powerMatrix");
		// System.out.println(Arrays.toString(powerMatrix[184]));
	}// Of revertPowerMatrix

	public void getModifiedMatrix() {
		modifiedMatrix = new double[numInstances()][numInstances()];
		double[][] EETNMatrix = new double[numInstances()][numInstances()];
		for (int i = 0; i < numInstances(); i++) {
			for (int j = 0; j < numInstances(); j++) {
				EETNMatrix[i][j] = (1 + 0.0) / numInstances();
			} // Of for j
		} // Of for i

		for (int i = 0; i < numInstances(); i++) {
			for (int j = 0; j < numInstances(); j++) {
				modifiedMatrix[i][j] = powerMatrix[i][j] * zuni + (1 - zuni) * EETNMatrix[i][j];
			} // Of for j
		} // Of for i
	}// Of getModifiedMatrix

	public static double getLevelValue(double[] vectorOne, double[] vectorTwo) {
		double sum = 0;
		for (int i = 0; i < vectorOne.length; i++) {
			sum += Math.abs(vectorOne[i] - vectorTwo[i]);
		} // Of for i

		return sum;
	}// Of getLevelValue

	public void getRankVector() {
		rankVector = new double[numInstances()];

		for (int i = 0; i < numInstances(); i++) {
			rankVector[i] = 1;
		} // Of for i

		times = 0;
		double currentOrderNorm = Double.MAX_VALUE;
		while (currentOrderNorm > threshold) {
			double[] tempRankVector = new double[numInstances()];
			for (int i = 0; i < numInstances(); i++) {
				for (int j = 0; j < numInstances(); j++) {
					tempRankVector[i] += modifiedMatrix[i][j] * rankVector[j];
				} // Of for j
			} // Of for i

			currentOrderNorm = 0;
			for (int i = 0; i < rankVector.length; i++) {
				currentOrderNorm += Math.abs(tempRankVector[i] - rankVector[i]);
			} // Of for i

			for (int i = 0; i < rankVector.length; i++) {
				rankVector[i] = tempRankVector[i];
			} // Of for i

			times++;
		} // Of while

		int[] indexRank = mergeSortToIndices(rankVector);

		System.out.println("排名序列号" + Arrays.toString(indexRank));

	}// Of getRankVector

	public void clusteringAndbuy() {

	}// clusteringAndbuy

	public RP[] getRepresentativeInstance() {
		RP[] reperesentativeSet = new RP[numInstances()];
		for (int i = 0; i < numInstances(); i++) {
			reperesentativeSet[i] = new RP();
			reperesentativeSet[i].instance = instance(i);
			reperesentativeSet[i].index = i;
			reperesentativeSet[i].rankValue = rankVector[i];
			reperesentativeSet[i].predictedLabel = -1;

		} // Of for i

		return reperesentativeSet;
	}// Of getRepresentativeInstance

	public int[][] computeBlockInformation(int tempBlocks, int[] tempclusterIndex) {

		// tempBlocks = tempclusterIndex.length;
		blockInformation = new int[tempBlocks][];

		for (int i = 0; i < tempBlocks; i++) {
			// Scan to see how many elements
			int tempElements = 0;
			for (int j = 0; j < tempclusterIndex.length; j++) {
				if (tempclusterIndex[j] == i) {
					tempElements++;
				} // Of if
			} // Of for k
				// Copy to the list
			blockInformation[i] = new int[tempElements];
			tempElements = 0;
			for (int j = 0; j < tempclusterIndex.length; j++) {
				if (tempclusterIndex[j] == i) {
					blockInformation[i][tempElements] = j;
					tempElements++;
				} // Of if
			} // Of for k
		} // Of for i
		System.out.println("information" + Arrays.deepToString(blockInformation));
		return blockInformation;
	}// Of computeBlockInforma

	public static RP[] getSortRepresentativeInstance(RP[] representativeSet) {
		RP tempRepresentative = new RP();
		for (int i = 0; i < representativeSet.length; i++) {
			for (int j = i + 1; j < representativeSet.length; j++) {
				if (representativeSet[i].rankValue < representativeSet[j].rankValue) {
					tempRepresentative = representativeSet[i];
					representativeSet[i] = representativeSet[j];
					representativeSet[j] = tempRepresentative;
				} // Of if
			} // Of for j
		} // Of for i

		return representativeSet;
	}// Of getSortRepresentativeInstance

	public int[] getCriticalInstances(double percentage) {

		// 重构；RP对象包括instance index 排名值
		RP[] representativeSet = getRepresentativeInstance();
		RP[] sortRepresentativeSet = getSortRepresentativeInstance(representativeSet);

		int count = (int) (percentage * sortRepresentativeSet.length);
		RP[] criticalSet = new RP[count];
		RP[] commonSet = new RP[sortRepresentativeSet.length - count];

		int countCritical = 0;
		int countCommon = 0;
		for (int i = 0; i < sortRepresentativeSet.length; i++) {
			if (i < count) {
				criticalSet[countCritical] = sortRepresentativeSet[i];
				countCritical++;
			} else {
				commonSet[countCommon] = sortRepresentativeSet[i];
				countCommon++;
			} // Of if
		} // Of for i

		// System.out.println("打印普通实例的标签");
		//
		// for (int i = 0; i < commonSet.length; i++) {
		//
		// System.out.println(""+commonSet[i].instance.classValue());
		// }

		// 重要实例之间的相似度
		// 在构树的时候是否可以可以按照距离来构建相似度？
		double[][] matrix = new double[criticalSet.length][criticalSet.length];
		for (int i = 0; i < criticalSet.length; i++) {
			Instance tempOne = criticalSet[i].instance;
			for (int j = 0; j < criticalSet.length; j++) {
				Instance tempTwo = criticalSet[j].instance;
				matrix[i][j] = similarityMatrix[criticalSet[i].index][criticalSet[j].index];
			} // Of for j
		} // Of for i

		System.out.println("重要实例的个数" + criticalSet.length);
		System.out.println("重要实例间的相似度" + Arrays.deepToString(matrix));

		// 将重要实例构建成一颗二叉树
		Tree binaryTree = new Tree(criticalSet.length);
		binaryTree.buildMatrix(matrix);
		// 建树
		binaryTree.buildTree();

		alreadyClassified = new boolean[criticalSet.length];
		// 买标签的个数
		numTeach = 0;
		// 预测标签的个数
		numPredict = 0;
		// 最有一步投票的个数
		numVote = 0;
		// 聚类阈值
		
		// 买和预测标签

		int timeWhile = 0;
		while (true) {
			// 遍历树、聚类
			timeWhile++;
			ArrayList<Integer> list = new ArrayList<Integer>();
			clusterNumbers = new int[criticalSet.length];
			for (int i = 0; i < clusterNumbers.length; i++) {
				clusterNumbers[i] = -1;
			}
			int tempClusterNum = 0;
			// 用于在递归中构建最大簇树
			binaryTree.clusterNum = 0;
			binaryTree.clusterNumbers = clusterNumbers;
			// 遍历左树
			binaryTree.clusteringTree(binaryTree.root.getLeftChild(), list, tempClusterNum, treethreshold);
			System.out.println("开始遍历右子树");
			// 遍历右树
			binaryTree.clusteringTree(binaryTree.root.getRightChild(), list, ++binaryTree.clusterNum, treethreshold);
			System.out.println("the clustering" + Arrays.toString(clusterNumbers));

			// 将聚类结果存入分块
			// 一共有多少类
			int clusterMax = getMax(clusterNumbers);
			computeBlockInformation(clusterMax + 1, clusterNumbers);

			boolean[] tempBlockProcessed = new boolean[clusterMax + 1];
			int tempUnProcessedBlocks = 0;
			for (int i = 0; i < blockInformation.length; i++) {
				tempBlockProcessed[i] = true;
				for (int j = 0; j < blockInformation[i].length; j++) {

					if (!alreadyClassified[blockInformation[i][j]]) {
						tempBlockProcessed[i] = false;
						tempUnProcessedBlocks++;
						break;
					} // of if
				} // of for j
			} // of for i

			for (int i = 0; i < blockInformation.length; i++) {
				// Step 2.3.1

				if (tempBlockProcessed[i]) {
					continue;
				} // of if

				if (blockInformation[i].length < 3) {

					for (int j = 0; j < blockInformation[i].length; j++) {
						if (!alreadyClassified[blockInformation[i][j]]) {
							if (numTeach >= 40) {
								break;
							} // of if

							criticalSet[blockInformation[i][j]].predictedLabel = (int) criticalSet[blockInformation[i][j]].instance
									.classValue();

							System.out.println("maidaode baoqianshi"
									+ (int) criticalSet[blockInformation[i][j]].instance.classValue());
							alreadyClassified[blockInformation[i][j]] = true;
							numTeach++;
							// System.out.println("numTeach first = " + numTeach);
						} // of if
					} // of for j
				} // of if

				int tempNumTeach = 0;
				for (int j = 0; j < blockInformation[i].length; j++) {
					if (alreadyClassified[blockInformation[i][j]]) {
						continue;
					} // of if
					if (numTeach >= 40) {
						break;
					} // of if
					criticalSet[blockInformation[i][j]].predictedLabel = (int) criticalSet[blockInformation[i][j]].instance
							.classValue();
					System.out.println(
							"maidaode baoqianshi" + (int) criticalSet[blockInformation[i][j]].instance.classValue());
					alreadyClassified[blockInformation[i][j]] = true;
					numTeach++;

					tempNumTeach++;

					if (tempNumTeach > 3) {
						break;
					} // of if

				} // of for j

			} // of for i
			boolean tempPure = true;

			for (int i = 0; i < blockInformation.length; i++) {

				if (tempBlockProcessed[i]) {
					continue;
				} // of if
				boolean tempFirstLable = true;

				int tempCurrentInstance;
				int tempLable = 0;

				for (int j = 0; j < blockInformation[i].length; j++) {
					tempCurrentInstance = blockInformation[i][j];
					if (alreadyClassified[tempCurrentInstance]) {

						if (tempFirstLable) {
							tempLable = criticalSet[blockInformation[i][j]].predictedLabel;
							tempFirstLable = false;
						} else {
							if (tempLable != criticalSet[blockInformation[i][j]].predictedLabel) {
								tempPure = false;
								break;
							} // of if
						} // of if
					} // of if
				} // of for j

				if (tempPure) {
					for (int j = 0; j < blockInformation[i].length; j++) {
						if (!alreadyClassified[blockInformation[i][j]]) {
							criticalSet[blockInformation[i][j]].predictedLabel = tempLable;
							alreadyClassified[blockInformation[i][j]] = true;
							numPredict++;
							// predictLabel[tempLable]++;
						} // of if
					} // of for j
				} // of if
			} // of for i

			treethreshold = treethreshold * 1.5;

			if (tempUnProcessedBlocks == 0) {
				break;
			} // of if

			if (numTeach >= 40) {
				break;
			} // of if

		} // of while
		System.out.println("timeWhile" + timeWhile);
		// 投票重要实例

		// 遍历树、聚类
		ArrayList<Integer> list = new ArrayList<Integer>();
		clusterNumbers = new int[criticalSet.length];

		for (int i = 0; i < clusterNumbers.length; i++) {
			clusterNumbers[i] = -1;
		}
		int tempClusterNum = 0;
		binaryTree.clusterNumbers = clusterNumbers;

		// 遍历左树
		binaryTree.clusteringTree(binaryTree.root.getLeftChild(), list, tempClusterNum, 0.37);
		System.out.println("开始遍历右子树");

		// 遍历右树
		binaryTree.clusteringTree(binaryTree.root.getRightChild(), list, ++binaryTree.clusterNum, 0.37);
		System.out.println("the clustering" + Arrays.toString(clusterNumbers));
		// 将聚类结果存入分块
		// 一共有多少类

		int clusterMax = getMax(clusterNumbers);
		computeBlockInformation(clusterMax + 1, clusterNumbers);

		int[][] vote = new int[clusterMax + 1][clusterMax + 1];
		int voteIndex = -1;

		for (int i = 0; i < blockInformation.length; i++) {

			for (int j = 0; j < blockInformation[i].length; j++) {

				for (int k = 0; k <= clusterMax; k++) {

					if (criticalSet[blockInformation[i][j]].predictedLabel == k) {
						vote[i][k]++;
					} // of if
				} // of for k
			} // of for j

			voteIndex = getMaxIndex(vote[i]);

			for (int j = 0; j < blockInformation[i].length; j++) {
				if (criticalSet[blockInformation[i][j]].predictedLabel == -1) {
					criticalSet[blockInformation[i][j]].predictedLabel = voteIndex;
					numVote++;
					// predictLabel[voteIndex]++;
				} // of if
			} // of for j

		} // of for i
		System.out.println("numVote = " + numVote + ",numTeach = " + numTeach);

		// 计算剩余实例与重要实例的相似度
		double[][] secondMatrix = new double[commonSet.length][criticalSet.length];

		for (int i = 0; i < commonSet.length; i++) {
			Instance tempOne = commonSet[i].instance;
			for (int j = 0; j < criticalSet.length; j++) {
				Instance tempTwo = criticalSet[j].instance;
				secondMatrix[i][j] = similarityMatrix[commonSet[i].index][criticalSet[j].index];
			} // Of for j
		} // Of for i

		// 需要计算重要实例的预测标签的最大值
		int maxClass = Integer.MIN_VALUE;
		;

		for (int i = 0; i < criticalSet.length; i++) {

			if (criticalSet[i].predictedLabel > maxClass)
				maxClass = criticalSet[i].predictedLabel;
		} // of for i
			// 计算剩余实例与簇之间的平均相似度
		for (int i = 0; i < commonSet.length; i++) {
			Instance instanceOne = commonSet[i].instance;
			double[] arrayAvergerSimilar = new double[maxClass + 1];
			for (int k = 0; k < maxClass + 1; k++) {
				int countOne = 0;
				for (int j = 0; j < criticalSet.length; j++) {
					if (criticalSet[j].predictedLabel == k) {
						Instance instanceTwo = criticalSet[j].instance;
						arrayAvergerSimilar[k] += getSimilarityValue(instanceOne, instanceTwo);
						countOne++;
					} // Of if
				} // of for j
				arrayAvergerSimilar[k] = (arrayAvergerSimilar[k] + 0.0) / countOne;

			} // Of for k
			int[] arraySortAvergerSimilar = new int[maxClass + 1];
			arraySortAvergerSimilar = mergeSortToIndices(arrayAvergerSimilar);
			commonSet[i].predictedLabel = arraySortAvergerSimilar[0];

		} // Of for i

		RP[] resultSet = new RP[numInstances()];
		for (int i = 0; i < count; i++) {
			resultSet[i] = criticalSet[i];
		} // Of for i
		for (int i = 0; i < numInstances() - count; i++) {
			resultSet[i + count] = commonSet[i];
		} // Of for i

		RP tempObject = null;
		for (int i = 0; i < numInstances(); i++) {
			for (int j = i; j < numInstances(); j++) {
				if (resultSet[i].index > resultSet[j].index) {
					tempObject = resultSet[i];
					resultSet[i] = resultSet[j];
					resultSet[j] = tempObject;
				} // Of if
			} // Of for j
		} // Of for i

		int[] predictedLabelVector = new int[numInstances()];
		for (int i = 0; i < numInstances(); i++) {
			predictedLabelVector[i] = resultSet[i].predictedLabel;
		} // Of for i

		return predictedLabelVector;
	}// Of getCriticalInstances

	public void delete(boolean[] paraIndices) {
		for (int i = numInstances() - 1; i >= 0; i--) {
			if (!paraIndices[i]) {
				delete(i);
			} // Of if
		} // Of for i
	}// Of delete

	public RDCnew[] divideInTwo(double paraPercentage) throws Exception {
		boolean[] firstInclusionArray = SimpleTool.generateBooleanArrayForDivision(numInstances(), paraPercentage);
		RDCnew firstDecisionSystem = new RDCnew(this);
		firstDecisionSystem.delete(firstInclusionArray);

		boolean[] secondInclusionArray = SimpleTool.revertBooleanArray(firstInclusionArray);
		RDCnew secondDecisionSystem = new RDCnew(this);
		secondDecisionSystem.delete(secondInclusionArray);

		RDCnew[] subsets = new RDCnew[2];
		subsets[0] = firstDecisionSystem;
		subsets[1] = secondDecisionSystem;

		return subsets;
	}// Of divideInTwo

	public int[] getRankDividedSets(double percentage) {
		// getSimilarityMatrix();

		getDistanceMatrix();
		getSimilarityMatrix();
		asignMostSimilar();
		// asignMostClose();
		revertPowerMatrix();
		getModifiedMatrix();
		getRankVector();
		int[] predictedLabelVector = this.getCriticalInstances(percentage);

		return predictedLabelVector;
	}// Of getRankDividedSets

	public RDCnew[] getRandomDividedSets(double percentage) {
		RDCnew[] randomDividedSets = null;
		try {
			randomDividedSets = divideInTwo(percentage);
		} catch (Exception e) {
			System.out.println(e.getStackTrace());
			System.out.println("Error occured in getRandomDividedSets(double).1");
		} // Of try

		return randomDividedSets;
	}// Of getRandomDividedSets

	public int[] getActualLabelVector() {
		int[] labelVector = new int[numInstances()];
		for (int i = 0; i < numInstances(); i++) {
			labelVector[i] = (int) instance(i).classValue();
		} // Of for i

		return labelVector;
	}// Of getActualLabelVector

	public static Instances generateUnsupervisedSet(RDCnew dataSet) {
		Instances initialSet = new Instances(dataSet);
		Instances clusterSet = new Instances(dataSet);
		clusterSet.delete();
		try {
			Remove filter = new Remove();
			filter.setAttributeIndices("" + (initialSet.classIndex() + 1));
			filter.setInputFormat(initialSet);
			clusterSet = weka.filters.Filter.useFilter(initialSet, filter);
		} catch (Exception e) {
			System.out.println(e.getStackTrace());
			System.out.println("Error occured in generateUnsupervisedSet().1");
		} // Of try

		return clusterSet;
	}// Of generateUnsupervisedSet

	public static void experiment() {
		RDCnew dataSet = null;
		try {
			FileReader fileReader = new FileReader(fileAddress);
			dataSet = new RDCnew(fileReader);
			fileReader.close();
			System.out.println("数据集大小" + dataSet.numInstances());
			dataSet.setClassIndex(dataSet.numAttributes() - 1);
		} catch (Exception e) {
			System.out.println("Error occured in experiment().1!");
			e.getStackTrace();
		} // Of try

		RDCnew[] tempSets = null;
		RDCnew usedSet = null;
		try {
			tempSets = dataSet.divideInTwo(0.99);
			// usedSet = tempSets[0];
			usedSet = dataSet;
		} catch (Exception e) {
			System.out.println("Error occured in experiment().2!");
			e.getStackTrace();
		} // OF try

		numClass = usedSet.attribute(dataSet.numAttributes() - 1).numValues();

		int[] actualLabelVector = usedSet.getActualLabelVector();
		System.out.println("真实标签" + Arrays.toString(actualLabelVector));
		System.out.println("chansfsu" + actualLabelVector.length);
		int[] rankPredictedLabelVector = usedSet.getRankDividedSets(percentage);
		System.out.println("预测标签" + Arrays.toString(rankPredictedLabelVector));
		getPredictionAccuracy(rankPredictedLabelVector, actualLabelVector);
		Instances clusterSet = generateUnsupervisedSet(usedSet);
		int[] KMLabelVector = simpleKMeansCluster(clusterSet);
		int[] EMLabelVector = EMCluster(clusterSet);
		int[] FFLabelVector = farthestFirstCluster(clusterSet);
		int[] HCLabelVector = hierachicalCluster(clusterSet);

		double[][] evaluationVector = new double[5][];

		evaluationVector[0] = getEvaluation(actualLabelVector, rankPredictedLabelVector);
		evaluationVector[1] = getEvaluation(actualLabelVector, KMLabelVector);
		evaluationVector[2] = getEvaluation(actualLabelVector, EMLabelVector);
		evaluationVector[3] = getEvaluation(actualLabelVector, FFLabelVector);
		evaluationVector[4] = getEvaluation(actualLabelVector, HCLabelVector);

		for (int i = 0; i < 5; i++) {
			System.out.println(Arrays.toString(evaluationVector[i]));
		} // Of for i

	}// Of experiment

	public static double[] getEvaluation(int[] actualVector, int[] predictedVector) {
		double[] evaluationVector = new double[3];
		int SS = 0;
		int SD = 0;
		int DS = 0;
		int DD = 0;
		for (int i = 0; i < actualVector.length; i++) {
			for (int j = i + 1; j < actualVector.length; j++) {
				if (actualVector[i] == actualVector[j]) {
					if (predictedVector[i] == predictedVector[j]) {
						SS++;
					} else {
						SD++;
					} // Of if
				} else {
					if (predictedVector[i] == predictedVector[j]) {
						DS++;
					} else {
						DD++;
					} // Of if
				} // Of if
			} // Of for j
		} // Of for i

		int length = actualVector.length;
		double JA = (SS + 0.0) / (SS + SD + DS);
		double FM = Math.sqrt(((SS + 0.0) / (SS + SD)) * ((SS + 0.0) / (SS + DS)));
		double RA = 2 * (SS + DD + 0.0) / (length * (length - 1));

		evaluationVector[0] = JA;
		evaluationVector[1] = FM;
		evaluationVector[2] = RA;

		return evaluationVector;
	}// Of getEvaluation

	public static void main(String[] args) {

		experiment();

	}// Of main

	public static int[] simpleKMeansCluster(Instances clusterSet) {
		int[] clusterLabel = new int[clusterSet.numInstances()];

		try {
			SimpleKMeans cluster = new SimpleKMeans();
			cluster.setNumClusters(numClass);
			cluster.buildClusterer(clusterSet);

			for (int i = 0; i < clusterSet.numInstances(); i++) {
				clusterLabel[i] = cluster.clusterInstance(clusterSet.instance(i));
			} // Of for i
		} catch (Exception e) {
			System.out.println("Something Error Occured in simpleKMeansCluster(Instances)!");
			e.getStackTrace();
		} // Of try

		return clusterLabel;
	}// Of simpleKMeansCluster

	public static int[] EMCluster(Instances clusterSet) {
		int[] clusterLabel = new int[clusterSet.numInstances()];

		try {
			String[] options = weka.core.Utils.splitOptions("-I 100 -N " + numClass + " -M 1.0E-6 -S 100");
			EM cluster = new EM();
			cluster.setOptions(options);
			cluster.buildClusterer(clusterSet);

			for (int i = 0; i < clusterSet.numInstances(); i++) {
				clusterLabel[i] = cluster.clusterInstance(clusterSet.instance(i));
			} // Of for i
		} catch (Exception e) {
			System.out.println("Something Error Occured in EMCluster(Instaces)!");
			e.getStackTrace();
		} // Of try

		return clusterLabel;
	}// Of EMCluster

	public static int[] farthestFirstCluster(Instances clusterSet) {
		int[] clusterLabel = new int[clusterSet.numInstances()];

		try {
			String[] options = new String[2];
			options[0] = "-S";
			options[1] = "100";
			FarthestFirst cluster = new FarthestFirst();
			cluster.setOptions(options);
			cluster.setNumClusters(numClass);
			cluster.buildClusterer(clusterSet);

			for (int i = 0; i < clusterSet.numInstances(); i++) {
				clusterLabel[i] = cluster.clusterInstance(clusterSet.instance(i));
			} // Of for i
		} catch (Exception e) {
			System.out.println("Something Error Occured in farthestFirstCluster(Instances)!");
			e.getStackTrace();
		} // Of try

		return clusterLabel;
	}// Of farthestFirstCluster

	public static int[] hierachicalCluster(Instances clusterSet) {
		int[] clusterLabel = new int[clusterSet.numInstances()];

		try {
			String[] options = new String[2];
			options[0] = "-L";
			options[1] = "WARD";
			HierarchicalClusterer cluster = new HierarchicalClusterer();
			cluster.setOptions(options);
			cluster.setNumClusters(numClass);
			cluster.buildClusterer(clusterSet);

			for (int i = 0; i < clusterSet.numInstances(); i++) {
				clusterLabel[i] = cluster.clusterInstance(clusterSet.instance(i));
			} // Of for i
		} catch (Exception e) {
			System.out.println("Something Error Occured in hierahicalCluster(Instances)!");
			e.getStackTrace();
		} // Of try

		return clusterLabel;
	}// Of hierachicalCluster

	public static int[] mergeSortToIndices(double[] paraArray) {
		int tempLength = paraArray.length;
		int[][] resultMatrix = new int[2][tempLength];// 两个维度交换存储排序tempIndex控制

		// Initialize
		int tempIndex = 0;
		for (int i = 0; i < tempLength; i++) {
			resultMatrix[tempIndex][i] = i;
		} // Of for i
			// System.out.println("Initialize, resultMatrix = " +
			// Arrays.deepToString(resultMatrix));

		// Merge
		int tempCurrentLength = 1;
		// The indices for current merged groups.
		// int tempFirstStart, tempSecondStart, tempSecondEnd;
		while (tempCurrentLength < tempLength) {
			// System.out.println("tempCurrentLength = " + tempCurrentLength);
			// Divide into a number of groups
			// Here the boundary is adaptive to array length not equal to 2^k.
			// ceil是向上取整函数
			// System.out.println("dddddddddddddMath.ceil(tempLength + 0.0 /
			// tempCurrentLength) / 2 "+Math.ceil((tempLength + 0.0) /
			// tempCurrentLength/2.0) +","+tempCurrentLength+"<"+tempLength);
			for (int i = 0; i < Math.ceil((tempLength + 0.0) / tempCurrentLength / 2.0); i++) {// 定位到哪一块

				// Boundaries of the group
				// System.out.println("for循环；第"+i+"块排序"+";当前每块长度tempCurrentLength
				// :"+tempCurrentLength);
				int tempFirstStart = i * tempCurrentLength * 2;
				// tempSecondStart定位第二块开始的位置index
				int tempSecondStart = tempFirstStart + tempCurrentLength;// 可以用于判断是否是最后一小块，并做初始化的工作
				int tempSecondEnd = tempSecondStart + tempCurrentLength - 1;
				if (tempSecondEnd >= tempLength) { // 控制最后一小块。若超过了整体长度，则当tempSecondEnd定位到数组最后
					tempSecondEnd = tempLength - 1;
				} // Of if
					// Merge this group
				int tempFirstIndex = tempFirstStart;
				int tempSecondIndex = tempSecondStart;
				int tempCurrentIndex = tempFirstStart;
				// System.out.println("Before merge");
				if (tempSecondStart >= tempLength) {
					for (int j = tempFirstIndex; j < tempLength; j++) {
						resultMatrix[(tempIndex + 1) % 2][tempCurrentIndex] = resultMatrix[tempIndex % 2][j];
						// tempFirstIndex++;
						tempCurrentIndex++;
					} // Of for j
					break;
				} // Of if

				while ((tempFirstIndex <= tempSecondStart - 1) && (tempSecondIndex <= tempSecondEnd)) {// 真正开始做排序的工作
					if (paraArray[resultMatrix[tempIndex % 2][tempFirstIndex]] >= paraArray[resultMatrix[tempIndex
							% 2][tempSecondIndex]]) {

						resultMatrix[(tempIndex + 1) % 2][tempCurrentIndex] = resultMatrix[tempIndex
								% 2][tempFirstIndex];

						tempFirstIndex++;
					} else {
						resultMatrix[(tempIndex + 1) % 2][tempCurrentIndex] = resultMatrix[tempIndex
								% 2][tempSecondIndex];

						tempSecondIndex++;
					} // Of if
					tempCurrentIndex++;

				} // Of while

				// Remaining part
				// System.out.println("Copying the remaining part");
				for (int j = tempFirstIndex; j < tempSecondStart; j++) {
					resultMatrix[(tempIndex + 1) % 2][tempCurrentIndex] = resultMatrix[tempIndex % 2][j];
					tempCurrentIndex++;

				} // Of for j
				for (int j = tempSecondIndex; j <= tempSecondEnd; j++) {
					resultMatrix[(tempIndex + 1) % 2][tempCurrentIndex] = resultMatrix[tempIndex % 2][j];
					tempCurrentIndex++;
				} // Of for j
					// paraArray=resultMatrix[0];
					// System.out.println("After copying remaining part");
					// System.out.println("Round " + tempIndex + ", resultMatrix = "
					// + Arrays.deepToString(resultMatrix));

			} // Of for i
				// System.out.println("Round " + tempIndex + ", resultMatrix = "
				// + Arrays.deepToString(resultMatrix));
			tempCurrentLength *= 2;
			tempIndex++;
		} // Of while
			// System.out.println("resultSortedIndices = " +
			// Arrays.toString(resultMatrix[tempIndex % 2]));

		return resultMatrix[tempIndex % 2];
	}// Of mergeSortToIndices

}// Of Class RDC

class RP {
	public Instance instance;
	public int index;
	public double rankValue;
	public int predictedLabel;
	public int type;
	public int master;
}// Of Class RP
