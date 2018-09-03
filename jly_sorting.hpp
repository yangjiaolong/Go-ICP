/********************************************************************
Sorting/Selection functions for the Go-ICP Algorithm
Last modified: Jan 27, 2015

"Go-ICP: Solving 3D Registration Efficiently and Globally Optimally"
Jiaolong Yang, Hongdong Li, Yunde Jia
International Conference on Computer Vision (ICCV), 2013

Copyright (C) 2013 Jiaolong Yang (BIT and ANU)

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
*********************************************************************/

#ifndef JLY_SORING_HPP
#define JLY_SORING_HPP

#define INTRO_K 5

#define INSERTION_NUM 5

//sorting in ascending order

template <typename T>
inline size_t median_of_st_mid_en(T * data, size_t st, size_t en)
{
	size_t mid = (st+en)/2;
	if(data[st] < data[en])
	{
		if(data[mid] < data[st])//median is data[st]
			return st;
		else if(data[mid] < data[en]) //median is data[md]
			return mid;
		else //median is data[en]
			return en;
	}
	else // data[en] <= data[st]
	{
		if(data[mid] < data[en])//median is data[en]
			return en;
		else if(data[mid] < data[st]) //median is data[md]
			return mid;
		else //median is data[st]
			return st;
	}
}

//median of 3 numbers
template <typename T>
inline size_t median_of_3(T * data, size_t st)
{
	T* data_ = data + st;

	if(data_[0] < data_[2])
	{
		if(data_[1] < data_[0])//median is data[0]
			return st;
		else if(data_[1] < data_[2]) //median is data[1]
			return st + 1;
		else //median is data[2]
			return st + 2;
	}
	else // data[2] <= data[0]
	{
		if(data[1] < data[2])//median is data[2]
			return st + 2;
		else if(data[1] < data[0]) //median is data[1]
			return st + 1;
		else //median is data[0]
			return st;
	}
}

//median of 5 numbers with 6 comparisons
template <typename T>
inline size_t median_of_5(T * data, size_t st)
{
	T* data_ = data + st;
	T tmp;

	if(data_[0] > data_[1])
	{
		tmp = data_[0];
		data_[0] = data_[1];
		data_[1] = tmp;
	}

	if(data_[2] > data_[3])
	{
		tmp = data_[2];
		data_[2] = data_[3];
		data_[3] = tmp;
	}

	if(data_[0] < data_[2])
	{
		tmp = data_[4];
		data_[4] = data_[0];

		if(tmp < data_[1])
			data_[0] = tmp;
		else
		{
			data_[0] = data_[1];
			data_[1] = tmp;
		}
	}
	else
	{
		tmp = data_[4];
		data_[4] = data_[2];

		if(tmp < data_[3])
			data_[2] = tmp;
		else
		{
			data_[2] = data_[3];
			data_[3] = tmp;
		}
	}

	if(data_[0] < data_[2])
	{
		if(data_[1] < data_[2])
			return st + 1;
		else
			return st + 2;
	}
	else
	{
		if(data_[0] < data_[3])
			return st;
		else
			return st + 3;
	}
}

template <typename T>
size_t median_of_medians(T * data, size_t st, size_t en)
{
	size_t l = en-st+1;
	size_t numof5 = l / 5;
	if(l % 5 != 0)
		numof5 ++;

	T tmp;
	size_t subst = st, suben = st + 4;
	size_t i, medind;

	//fist (numof5 - 1) groups
	for(i = 0; i < numof5 - 1; i++, subst += 5, suben += 5)
	{
		medind = median_of_5(data, subst);

		tmp = data[st+i];
		data[st+i] = data[medind];
		data[medind] = tmp;
	}

	//last group
	{
		switch(en-subst+1)
		{
		case 3: // 3 elements
		case 4: // 4 elements
			medind = median_of_3(data, subst);
			break;
		case 5: // 5 elements
			medind = median_of_5(data, subst);
			break;
		default: // 1 or 2 elements
			medind = subst;
			break;
		}

		tmp = data[st+i];
		data[st+i] = data[medind];
		data[medind] = tmp;
	}
	
	//median of medians
	if(numof5 > 5)
		return median_of_medians(data, st, st + numof5-1);
	else
	{
		switch(numof5)
		{
		case 3: // 3 elements
		case 4: // 4 elements
			return median_of_3(data, st);
			break;
		case 5: // 5 elements
			return median_of_5(data, st);
			break;
		default: // 1 or 2 elements
			return st;
			break;
		}
	}
}

template <typename T>
void insertion_sort(T * data, size_t st, size_t en)
{
	T tmp;
	size_t i, j;
	for(i = st+1; i <= en; i++)
		for(j = i; j > st && data[j-1] > data[j]; j--)
		{
			tmp = data[j-1];
			data[j-1] = data[j];
			data[j] = tmp;
		}
}

// Sort the given array in ascending order
// Stop immediately after the array is splitted into k small numbers and n-k large numbers
template <typename T>
void intro_select(T * data, size_t st, size_t en, size_t k)
{	
	T pivot;
	T tmp;

	//for(; st < en && data[st] > 0; st++);

	size_t l_pre = en-st+1;
	size_t l;
	size_t medind;

	bool quickselect = true;

	size_t i = 0;
	while(1)
	{
		if(st >= en)
			break;

		if(en - st <= INSERTION_NUM)
		{	
			insertion_sort(data, st, en);
			return;
		}

		if(quickselect && i++ == INTRO_K)
		{
			// switch to 'median of medians' if INTRO_K partations of quickselect fail to half the size 
			l = en-st+1;
			if(l*2 > l_pre)
				quickselect = false;
			
			l_pre = l;
			i = 0;
		}

		if(quickselect)
			//medind = st;
			medind = median_of_st_mid_en(data, st, en);
		else
			medind = median_of_medians(data, st, en);

		if(medind != st)
		{
			tmp = data[st];
			data[st] = data[medind];
			data[medind] = tmp;
		}

		size_t p = st;
		size_t left = st+1;
		size_t right = en;
		pivot = data[p];

		while(1)
		{
			while (left < right && pivot >= data[left])
				++left;
			while (left < right && pivot <= data[right])
				--right;

			if (left >= right)
				break;

			//swap left & right
			tmp = data[left];
			data[left] = data[right];
			data[right] = tmp;
		}

		size_t s = left-1;
		if(data[left] < pivot)
			s = left;
		//swap p & s
		data[p] = data[s];
		data[s] = pivot;

		if(s < k)
			st = s+1;
		else if(s > k)
			en = s-1;
		else //s == k
			break;
	}
}

#endif
