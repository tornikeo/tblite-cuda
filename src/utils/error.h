#ifndef UTILS_ERROR_H
#define UTILS_ERROR_H

/* 
   !> Error message
   type :: error_type

      !> Error code
      integer :: stat

      !> Payload of the error
      character(len=:), allocatable :: message

   end type error_type */

class error_type
{
  public:
  int stat {0};
  char message[32] {0};
};

#endif