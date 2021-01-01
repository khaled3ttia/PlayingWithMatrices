
int binarySearch(std::vector<int> &searchArray, int start, int end, int value){
	
	while (start <= end){
		
					
		int midPoint = start + (int)ceil((end - start) / 2);

		if (searchArray[midPoint] == value){
			return midPoint;
		}
	
		if (value > searchArray[midPoint]){

			start = midPoint + 1;

		}else {
			end = midPoint - 1;
		}
		

	}
	return -1;


}
