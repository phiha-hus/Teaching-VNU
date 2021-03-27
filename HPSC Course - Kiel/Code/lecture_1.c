/* 
Lecture 1: Begin with C, 
learn for, while, if ... else, switch
learn IO
learn data type: int, double, float, char 
*/
# include<stdio.h>

int main()

{
    int i;
    for (i=0 ; i < 20 ; i++);
        {
            printf("Gia tri cua so %d \n",i);
            if (5<i && i<12)
                {
                    printf("Thay Phi dep trai qua di mat ! \n");
                }
            else
            {
                printf("Em yeu thay Phi nhat qua dat \n");
            }
        }           
    
    printf("For loop ends here!!!!!!!!") ;
   
    i = 5 ;
    printf(i) ;
    do  {
            if (i>7)
            {
                printf("So i la: %d \n",i);
                i++;
            }
        }while (i< 10) ;

    float j;
    printf("Input a real number: \n");
    scanf("%f",&j);
    printf("So j can tim chinh la %f \n",j);

    return 0;
}