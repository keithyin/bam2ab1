
pub fn binary_search_lower_bound<T: Ord>(arr: &[T], target: &T) -> usize {
    let mut left = 0;
    let mut right = arr.len();
    while left < right {
        let mid = left + (right - left) / 2;
        if &arr[mid] < target {
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    left
}