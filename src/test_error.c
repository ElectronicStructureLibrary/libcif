/*
 Copyright (C) 2015 Y. Pouillon

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or 
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/

/**
 * @file test_error.c
 * @brief checks error handler
 */

#include <stdlib.h>
#include <stdio.h>

#include "vector3.h"
#include "test_common.h"


int main(void) {
  bool eid = 0;
  Vector3 v1, v2;

  /* Display basic information */
  DEBUG_PRINT("%s - test_error\nReport bugs to %s\n\n", PACKAGE_STRING,
    PACKAGE_BUGREPORT);
  DEBUG_PRINT("=== BEGIN test_error ===\n\n");

  /* FIXME: dummy check */
  DEBUG_PRINT("test_error: checking vector equality\n");
  DEBUG_PRINT("test_error: at the beginning, status = %d\n", eid);
  v3_new(&v1, 1.0, 2.0, 3.0);
  v3_new(&v2, 3.0, 2.0, 1.0);
  eid = v3_eq(&v1, &v2, 0.5);
  DEBUG_PRINT("test_error: after v3_eq, status = %d\n", eid);
  DEBUG_PRINT("\n");

  DEBUG_PRINT("=== END test_error ===\n");

  return 0;
}
