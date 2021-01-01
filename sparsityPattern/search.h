
int binarySearch(std::vector<int> &searchArray, int start, int end, int value){
	
	while (start <= end){
		
					
		int midPoint = start + (int)ceil((end - start) / 2);
		
		// The use for this algorithm in our tool is only to search within index arrays
		// in other parts of the application, we delete entries from index array by assigning
		// their values to -1 (since -1 is an invalid value for an array index)
		// So, we need to check first against the value of -1 

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
