#include <stdio.h>
main()
{
    int c;
    int inQuote = 0;
    while((c=getchar())!=EOF)
    {
	if(c=='"') inQuote = 1-inQuote;
	if(inQuote && c==',') putchar(' ');
	else putchar(c);
    }
}
