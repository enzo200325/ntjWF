/**
 * Author: 
 * Date: 
 * License: 
 * Source: 
 * Description: Set with operations to find the element by index \- find\_by\_order(i) \- and to find the index of an element \- \- 
 * Time: All operations in O(\log N).
 * Usage: ordered_set<int> seta; 
 * Status: should be good to go
 */

#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace __gnu_pbds;

template <typename T>
using ordered_set =
    tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;

template <typename T, typename U>
using ordered_map = tree<T, U, less<T>, rb_tree_tag, tree_order_statistics_node_update>;
